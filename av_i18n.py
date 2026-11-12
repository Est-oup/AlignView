"""French/English presentation; biological data and command arguments stay intact."""
from __future__ import annotations

from contextvars import ContextVar
import json
import os
import re
import xml.etree.ElementTree as ET

from av_locale_en import EN

_LANGUAGE=ContextVar('alignview_language',default='en' if os.environ.get('ALIGNVIEW_LANG')=='en' else 'fr')
_FOLDED={key.casefold():value for key,value in EN.items()}
_PATTERNS=[]
for key,value in EN.items():
    if '{' not in key:continue
    chunks=re.split(r'(\{[A-Za-z_][A-Za-z_0-9]*\})',key)
    names=[];parts=[]
    for chunk in chunks:
        if re.fullmatch(r'\{[A-Za-z_][A-Za-z_0-9]*\}',chunk):
            names.append(chunk[1:-1]);parts.append('(.+?)')
        else:parts.append(re.escape(chunk))
    _PATTERNS.append((re.compile('^'+''.join(parts)+'$',re.IGNORECASE),names,value))


def language():return _LANGUAGE.get()


def set_language(value):
    if value not in ('fr','en'):raise ValueError('Language must be fr or en')
    _LANGUAGE.set(value)


def translate(text,_depth=0):
    text=str(text)
    if language()!='en' or not text or _depth>8:return text
    if '\n' in text:return '\n'.join(translate(line,_depth+1) for line in text.split('\n'))
    core=text.strip()
    if not core:return text
    result=EN.get(core) or _FOLDED.get(core.casefold()) or _FOLDED.get(' '.join(core.split()).casefold())
    if result is None:
        if core.startswith('. '):result='. '+translate(core[2:],_depth+1)
    if result is None:
        for pattern,names,target in _PATTERNS:
            found=pattern.fullmatch(core)
            if found:
                result=target.format(**dict(zip(names,found.groups())))
                break
    if result is None:
        for separator in (' › ',' · ',': '):
            if separator in core:
                result=separator.join(translate(part,_depth+1) for part in core.split(separator))
                break
    if result is None:return text
    if core.isupper():result=result.upper()
    return text[:len(text)-len(text.lstrip())]+result+text[len(text.rstrip()):]


def localize_svg(svg):
    if language()!='en':return svg
    root=ET.fromstring(svg)
    for node in root.iter():
        if set(node.get('class','').split()) & {'tree-tip','tree-support','data-label'}:continue
        if node.tag.rsplit('}',1)[-1] in ('text','title','desc') and node.text:
            node.text=translate(node.text)
        if 'aria-label' in node.attrib:node.set('aria-label',translate(node.get('aria-label')))
    ET.register_namespace('','http://www.w3.org/2000/svg')
    root.set('lang','en')
    return ET.tostring(root,encoding='unicode')


def figure_title(svg,fallback='Figure'):
    root=ET.fromstring(svg)
    title=root.find('.//{http://www.w3.org/2000/svg}text')
    return ''.join(title.itertext()) if title is not None else translate(fallback)


def localize_html(page):
    if language()!='en':return page
    catalog=json.dumps(EN,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    script=r'''<script>(function(){
const words=__CATALOG__,folded=Object.fromEntries(Object.entries(words).map(([a,b])=>[a.toLowerCase(),b]));
const escape=x=>x.replace(/[.*+?^${}()|[\]\\]/g,'\\$&');
const patterns=Object.entries(words).filter(([a])=>a.includes('{')).map(([source,target])=>{
  const names=[];const expression=source.split(/(\{[A-Za-z_][A-Za-z_0-9]*\})/).map(part=>{
    if(/^\{[A-Za-z_][A-Za-z_0-9]*\}$/.test(part)){names.push(part.slice(1,-1));return '(.+?)';}return escape(part);
  }).join('');return [new RegExp('^'+expression+'$','i'),names,target];
});
function tr(text,depth=0){
 if(!text||depth>8)return text;
 if(text.includes('\n'))return text.split('\n').map(line=>tr(line,depth+1)).join('\n');
 const core=text.trim();if(!core)return text;let result=words[core]||folded[core.toLowerCase()];
 if(!result)for(const [pattern,names,target] of patterns){const found=core.match(pattern);if(found){
   const values=Object.fromEntries(names.map((name,i)=>[name,found[i+1]]));result=target.replace(/\{([A-Za-z_][A-Za-z_0-9]*)\}/g,(_,name)=>values[name]);break;
 }}
 if(!result)for(const separator of [' › ',' · ',': '])if(core.includes(separator)){result=core.split(separator).map(x=>tr(x,depth+1)).join(separator);break;}
 if(!result&&core.startsWith('. '))result='. '+tr(core.slice(2),depth+1);
 if(!result)return text;return text.slice(0,text.indexOf(core))+result+text.slice(text.indexOf(core)+core.length);
}
function update(){
 const walk=document.createTreeWalker(document.body,NodeFilter.SHOW_TEXT);const nodes=[];
 while(walk.nextNode())nodes.push(walk.currentNode);
 for(const node of nodes){
   if(node.parentElement?.closest('script,style,pre,.tree-tip,.tree-support,.data-label,.msa-name,.msa-seq,.tablewrap td:first-child,#vsearch td:nth-child(2)'))continue;
   const next=tr(node.nodeValue);if(next!==node.nodeValue)node.nodeValue=next;
 }
 for(const node of document.querySelectorAll('[placeholder],[aria-label],button[title]'))for(const attr of ['placeholder','aria-label','title']){
   if(node.hasAttribute(attr)){const old=node.getAttribute(attr),next=tr(old);if(old!==next)node.setAttribute(attr,next);}
 }
}
document.addEventListener('DOMContentLoaded',()=>{update();new MutationObserver(update).observe(document.body,{childList:true,subtree:true,characterData:true});});
})();</script>'''.replace('__CATALOG__',catalog)
    return page.replace('<html lang="fr">','<html lang="en">',1).replace('</body>',script+'</body>',1)
