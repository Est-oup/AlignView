"""Read FastTree Newick and draw a standalone SVG using the standard library."""
from __future__ import annotations

from dataclasses import dataclass, field
import html
import math


TREE_FIGURE = '45_fasttree_tree'


@dataclass(eq=False)
class TreeNode:
    label: str = ''
    length: float | None = None
    children: list[TreeNode] = field(default_factory=list)


def _tokens(source):
    """Keep quoted labels intact; ignore whitespace and Newick comments."""
    i = 0
    while i < len(source):
        char = source[i]
        if char.isspace():
            i += 1
        elif char == '[':
            depth = 1
            i += 1
            while i < len(source) and depth:
                if source[i] == '[': depth += 1
                elif source[i] == ']': depth -= 1
                i += 1
            if depth: raise ValueError('Invalid Newick: unterminated comment.')
        elif char == "'":
            label = []
            i += 1
            while i < len(source):
                if source[i] == "'":
                    i += 1
                    if i == len(source) or source[i] != "'": break
                label.append(source[i])
                i += 1
            else:
                raise ValueError('Invalid Newick: unterminated quoted label.')
            yield 'label', ''.join(label)
        elif char in '(),:;':
            yield char, char
            i += 1
        else:
            start = i
            while i < len(source) and not source[i].isspace() and source[i] not in "[]'(),:;":
                i += 1
            if i == start: raise ValueError('Invalid Newick: unexpected character.')
            yield 'label', source[start:i]


def parse_newick(source: str) -> TreeNode:
    """Parse one tree without recursion, including multifurcations and supports."""
    tokens = list(_tokens(source))
    position = 0

    def suffix(node, internal=False):
        nonlocal position
        if internal and position < len(tokens) and tokens[position][0] == 'label':
            node.label = tokens[position][1]
            position += 1
        if position < len(tokens) and tokens[position][0] == ':':
            position += 1
            if position == len(tokens) or tokens[position][0] != 'label':
                raise ValueError('Invalid Newick: missing branch length.')
            try:
                node.length = float(tokens[position][1])
            except ValueError as exc:
                raise ValueError('Invalid Newick: non-numeric branch length.') from exc
            if not math.isfinite(node.length) or node.length < 0:
                raise ValueError('Invalid Newick: negative or non-finite branch length.')
            position += 1

    root = None
    parents = []
    expect_node = True
    while position < len(tokens):
        kind, value = tokens[position]
        position += 1
        if expect_node:
            if kind not in ('(', 'label') or (kind == 'label' and not value):
                raise ValueError('Invalid Newick: expected a sequence or subtree.')
            node = TreeNode(label=value if kind == 'label' else '')
            if parents: parents[-1].children.append(node)
            elif root is None: root = node
            else: raise ValueError('Invalid Newick: plusieurs racines.')
            if kind == '(':
                parents.append(node)
            else:
                suffix(node)
                expect_node = False
        elif kind == ',' and parents:
            expect_node = True
        elif kind == ')' and parents:
            suffix(parents.pop(), internal=True)
        elif kind == ';' and not parents and position == len(tokens):
            return root
        else:
            raise ValueError('Invalid Newick: unexpected separator or parenthesis.')
    raise ValueError('Invalid Newick: incomplete tree (parentheses or semicolon).')


def tree_text(source: str) -> str:
    """Readable topology; numerical lengths are labels, not terminal spacing."""
    root=parse_newick(source)
    lines=['Tree topology (spacing is not to scale).',
           'Lengths shown: substitutions/site. Arbitrary drawing root.',
           'FastTree local support: not bootstraps. SVG figure: Save and share.', '']
    stack=[(root,'',True,True)]
    while stack:
        node,prefix,last,is_root=stack.pop()
        label=node.label or ('●' if node.children else '?')
        if node.length is not None:label+=': '+format(node.length,'.6g')
        lines.append(prefix+('' if is_root else '└─ ' if last else '├─ ')+label)
        child_prefix=prefix+('' if is_root else '   ' if last else '│  ')
        for index in range(len(node.children)-1,-1,-1):
            stack.append((node.children[index],child_prefix,index==len(node.children)-1,False))
    return '\n'.join(lines)


def tree_svg(source: str) -> str:
    """Rectangular phylogram; absent/all-zero lengths use an explicit cladogram."""
    root = parse_newick(source)
    nodes = []
    stack = [(root, None, 0, 0.0)]
    while stack:
        node, parent, depth, distance = stack.pop()
        if not math.isfinite(distance):
            raise ValueError('Invalid Newick: non-finite sum of branch lengths.')
        nodes.append((node, parent, depth, distance))
        stack.extend((child, node, depth + 1, distance + (child.length or 0.0))
                     for child in reversed(node.children))
    leaves = [node for node, _, _, _ in nodes if not node.children]
    if len(leaves) < 2:
        raise ValueError('A tree must contain at least two sequences.')
    max_distance = max(distance for _, _, _, distance in nodes)
    scaled = all(node.length is not None for node, parent, _, _ in nodes if parent) and max_distance > 0
    maximum = max_distance if scaled else max(depth for _, _, depth, _ in nodes)
    left, top, plot_width, row = 42, 114, 720, 30
    label_width = max(260, max(len(node.label) for node in leaves) * 8 + 30)
    width = left + plot_width + label_width + 30
    height = top + (len(leaves) - 1) * row + 116
    xs = {node: left + ((distance / maximum) if scaled else (depth / maximum)) * plot_width
          for node, _, depth, distance in nodes}
    ys = {node: top + i * row for i, node in enumerate(leaves)}
    for node, _, _, _ in reversed(nodes):
        if node.children:
            ys[node] = (ys[node.children[0]] + ys[node.children[-1]]) / 2
    esc = html.escape
    layout = 'phylogram' if scaled else 'cladogram'
    note = ('Horizontal lengths: substitutions per site.' if scaled else
            'Cladogram: missing or all-zero branch lengths; distances are not represented.')
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'role="img" aria-label="Arbre exploratoire FastTree" data-layout="{layout}" data-leaves="{len(leaves)}">',
        '<title>Exploratory tree · FastTree</title>',
        '<desc>Branchs and labels from the Newick file. Arbitrary drawing root, '
        'without biological rooting. FastTree support is local (SH-like), not bootstrap support.</desc>',
        f'<rect width="{width}" height="{height}" rx="14" fill="#102039"/>',
        '<g font-family="Arial,sans-serif">',
        '<text x="28" y="35" fill="#dce9f8" font-size="21" font-weight="bold">Exploratory tree · FastTree</text>',
        f'<text x="28" y="60" fill="#9fb4ce" font-size="12">{len(leaves)} sequences · Arbitrary drawing root, without biological rooting</text>',
        f'<text x="28" y="81" fill="#9fb4ce" font-size="12">{note}</text>',
    ]
    for node, parent, _, _ in nodes:
        x, y = xs[node], ys[node]
        if node.children:
            parts.append(f'<path d="M{x:.3f} {ys[node.children[0]]:.3f}V{ys[node.children[-1]]:.3f}" '
                         'stroke="#75a9ec" stroke-width="1.7" fill="none"/>')
        if parent:
            length = '' if node.length is None else f'{node.length:.12g}'
            parts.append(f'<path class="tree-branch" data-length="{length}" '
                         f'd="M{xs[parent]:.3f} {y:.3f}H{x:.3f}" stroke="#75a9ec" '
                         f'stroke-width="1.7" fill="none"><title>Branch: {length or "not specified"}</title></path>')
        if not node.children:
            parts.extend([
                f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" fill="#5bdbbb"/>',
                f'<text class="tree-tip" x="{x+10:.3f}" y="{y+4:.3f}" font-size="13" fill="#deebf9">{esc(node.label)}</text>',
            ])
        elif node.label:
            # A backed label on the incoming branch remains readable even when
            # a descendant has an almost-zero length and sits next to the node.
            label_width = len(node.label) * 6.5 + 12
            center = (xs[parent] + x) / 2 if parent else x + label_width / 2 + 6
            parts.append(f'<rect x="{center-label_width/2:.3f}" y="{y-9:.3f}" '
                         f'width="{label_width}" height="18" rx="4" fill="#102039"/>')
            parts.append(f'<text class="tree-support" x="{center:.3f}" y="{y+4:.3f}" '
                         f'text-anchor="middle" font-size="11" fill="#ffca7c">{esc(node.label)}</text>')
    footer = height - 70
    if scaled:
        # A 1/2/5 scale no longer than one fifth of the plot; normalized first
        # so even very small branch lengths do not overflow a pixel multiplier.
        target = maximum / 5
        if target > 0:
            magnitude = max(math.ulp(0.0), 10 ** math.floor(math.log10(target)))
            scale = next(value * magnitude for value in (5, 2, 1) if value * magnitude <= target)
            end = left + (scale / maximum) * plot_width
            parts.append(f'<path class="tree-scale" data-value="{scale:.12g}" d="M{left} {footer-5}v10m0 -5H{end:.3f}m0 -5v10" '
                         'stroke="#dce9f8" stroke-width="1.5" fill="none"/>')
            parts.append(f'<text x="{end+12:.3f}" y="{footer+4}" fill="#c2d4e8" font-size="12">{scale:.3g} substitutions/site</text>')
    parts.extend([
        f'<text x="28" y="{height-35}" fill="#9fb4ce" font-size="12">Internal labels: local FastTree support (SH-like), not bootstrap values.</text>',
        f'<text x="28" y="{height-15}" fill="#9fb4ce" font-size="11">Exploratory interpretation; depends on MSA quality and the selected model.</text>',
        '</g></svg>',
    ])

    return str('\n'.join(parts))
