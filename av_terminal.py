"""Small terminal primitives shared by menus, the pager and the MSA viewer."""
from __future__ import annotations

from contextlib import contextmanager
import os
import re
import select
import shutil
import sys
import termios
import tty
import unicodedata

SGR = re.compile(r'\x1b\[[0-9;]*m')
_screens = []


def fit(text, width):
    """Clip by terminal cells, preserving only SGR styling; never wrap a frame."""
    result = []
    cells = 0
    pos = 0
    text = str(text)
    while pos < len(text):
        match = SGR.match(text, pos)
        if match:
            result.append(match.group())
            pos = match.end()
            continue
        char = text[pos]
        pos += 1
        if unicodedata.category(char).startswith('C'):
            continue
        size = 0 if unicodedata.combining(char) else 2 if unicodedata.east_asian_width(char) in ('W', 'F') else 1
        if cells + size > max(0, width):
            break
        result.append(char)
        cells += size
    return ''.join(result) + ('\x1b[0m' if '\x1b[' in text else '')


def draw(lines):
    from av_i18n import translate
    lines=[translate(line) for line in lines]
    size = shutil.get_terminal_size((104, 32))
    frame = '\x1b[K\r\n'.join(fit(line, size.columns - 1) for line in lines[:max(1, size.lines - 1)])
    sys.stdout.write('\x1b[H' + frame + '\x1b[K\x1b[J')
    sys.stdout.flush()


@contextmanager
def screen(mouse=False):
    fd = sys.stdin.fileno()
    original = termios.tcgetattr(fd)
    try:
        # TCSAFLUSH discards fast input typed during a menu transition.
        tty.setcbreak(fd,when=termios.TCSANOW)
        prefix='\x1b[?1049h\x1b[?25l' if not _screens else ''
        if mouse and not any(_screens):prefix+='\x1b[?1000h\x1b[?1006h'
        _screens.append(mouse)
        sys.stdout.write(prefix)
        sys.stdout.flush()
        yield original
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, original)
        _screens.pop()
        suffix='\x1b[?1000l\x1b[?1006l' if mouse and not any(_screens) else ''
        suffix+='\x1b[0m\x1b[?25h\x1b[?1049l' if not _screens else '\x1b[?25l'
        sys.stdout.write(suffix)
        sys.stdout.flush()


def read_key(timeout=0.3):
    """Read a complete UTF-8 character or CSI/SS3 sequence, including modifiers."""
    fd = sys.stdin.fileno()
    if not select.select([fd], [], [], timeout)[0]:
        return 'tick'
    raw = os.read(fd, 1)
    if not raw:
        return 'eof'
    if raw == b'\x03':
        raise KeyboardInterrupt
    if raw != b'\x1b':
        size = 1 if raw[0] < 128 else 2 if raw[0] < 224 else 3 if raw[0] < 240 else 4
        for _ in range(size - 1):
            if not select.select([fd], [], [], .1)[0]:
                break
            raw += os.read(fd, 1)
        return raw.decode('utf-8', 'replace')
    if not select.select([fd], [], [], .08)[0]:
        return 'esc'
    prefix = os.read(fd, 1)
    if prefix not in (b'[', b'O'):
        return 'unknown'
    tail = b''
    for _ in range(40):
        if not select.select([fd], [], [], .08)[0]:
            break
        byte = os.read(fd, 1)
        tail += byte
        if byte and 0x40 <= byte[0] <= 0x7e:
            break
    code = tail.decode('ascii', 'replace')
    if code.startswith('<'):
        return 'mouse:' + code[1:]
    direct = {'A':'up', 'B':'down', 'C':'right', 'D':'left', 'H':'home', 'F':'end',
              '5~':'pageup', '6~':'pagedown', '1~':'home', '4~':'end', '7~':'home', '8~':'end',
              'P':'F1', '11~':'F1', 'Q':'F2', '12~':'F2'}
    if code in direct:
        return direct[code]
    match = re.fullmatch(r'(?:1|5|6);([2-8])([ABCDHF~])', code)
    if match:
        modifier, end = match.groups()
        direction = direct.get(end, 'pageup' if code.startswith('5;') else 'pagedown')
        return ('shift-' if modifier == '2' else 'ctrl-') + direction
    return 'unknown'
