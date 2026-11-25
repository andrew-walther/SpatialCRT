"""Inline the CTJ exhibit macros at their call sites so latexdiff can mark captions and cells."""
import re,sys
man,exh,out=sys.argv[1:4]
m=open(man).read(); e=open(exh).read()
defs={}; keep=[]; i=0
pat=re.compile(r'\\newcommand\{\\(\w+)\}(\[\d\])?\{')
pos=0; rest=[]
while True:
    mm=pat.search(e,pos)
    if not mm: rest.append(e[pos:]); break
    rest.append(e[pos:mm.start()])
    j=mm.end(); depth=1
    while depth:
        c=e[j]
        if c=='\\': j+=2; continue
        if c=='{': depth+=1
        elif c=='}': depth-=1
        j+=1
    name=mm.group(1); body=e[mm.end():j-1]
    if mm.group(2): keep.append(e[mm.start():j])  # parameterized helpers (CTJCat, CTJThumb) stay as macros
    else: defs[name]=body.lstrip('%\n')
    pos=j
helpers='\n'.join(keep)+'\n'
assert m.count('\\input{CTJ_Exhibits.tex}')==1
m=m.replace('\\input{CTJ_Exhibits.tex}',helpers)
for name,body in defs.items():
    m=re.sub(r'\\'+name+r'(?![A-Za-z])',lambda _:body,m)
open(out,'w').write(m)
print('inlined',sorted(defs),'kept',len(keep))

# Also inline the abstract macro defined in the manuscript preamble.
m=open(out).read()
k=m.index('\\newcommand{\\CTJAbstract}{%'); j=k+len('\\newcommand{\\CTJAbstract}{%'); d=1; p=j
while d:
    c=m[p]
    if c=='\\': p+=2; continue
    d+= (c=='{') - (c=='}'); p+=1
body=m[j:p-1]
m=m[:k]+m[p:]
n=len(re.findall(r'\\CTJAbstract(?![A-Za-z])',m))
m=re.sub(r'\\CTJAbstract(?![A-Za-z])',lambda _:body,m)
open(out,'w').write(m); print('abstract calls inlined',n)
