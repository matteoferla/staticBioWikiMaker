#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Written for python 3, not tested under 2.
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "13 / 12 / 15"
__doc__='''
staticBioWikiMaker(markdown, proteomics,entrytemplate,categorytemplate, indexpage) makes a static set of webpages based on a markdown file.

The level 1 headers on the markdown file are splint into entries and inserted into a template (entrytemplate), which has a $title, $menu and $content keywords for substitution.

The level 2 headers are parse differently. Members should contain gene symbols prefixed with a dollar sign. The script will extract the proteomics data from the proteomics file and median normalise it and show on that entry the genes tagged in that manner.
Tags contains keywords with a dollar sigil that are common to other entries. The site will generate a navigation box for each keyword with links to all entries which mention it.

The index page is an SVG with a layer called "clickareas" in front of the main picture and containing polygons of opacity 0% with ids matching the entries.
The script will link them and change the opacity onhover.
The entries will generate the dropdown menu $menu, but entries can be generated that are hardcoded in the template, such as an about page.

There are a few extra features not discussed.

Also to get google drive to serve the html as such the <id> needs to be put in https://googledrive.com/host/<id>
'''


import markdown as md
import re
import json
from string import Template
from collections import defaultdict
import xml.etree.ElementTree as et
import pandas as pd

from collections import OrderedDict, Callable

class DefaultOrderedDict(OrderedDict):  #http://stackoverflow.com/questions/6190331/can-i-do-an-ordered-default-dict-in-python
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,OrderedDict.__repr__(self))


class mdentry():
    plan_options=("increase", "decrease", "heterologous", "knockout", "mutate", "feed")  #sloppy enum
    def __init__(self, txt):
        self.title=re.match("\# (.*?)\n",txt).groups(1)[0]
        print(self.title)
        self.tags=[]
        self.compound=[]
        self.genes=[]
        self.summary="Summary not available"
        plandex=defaultdict(list)
        for s in txt.split("## "):
            subheader=re.match("\s*(.*?)\n",s).groups(1)[0]  #check that + is zero or more
            for element in "tags compounds genes".split():
                if subheader.lower() == element:
                    setattr(self,element,re.findall("\$(\w+)",s))
            if subheader.lower() == "summary":
                    self.summary=s.replace("Summary\n","")
            if subheader.lower() == "strategies":
                for plan in re.findall("\$(\w+)",s):
                    for typology in self.plan_options:
                        rex=re.match(typology+"_(\w+)",plan)
                        if rex:
                            plandex[typology].append(rex.groups(1)[0])
                txt=txt.replace(s,"Strategies\n"+" ".join(["$"+x for x in plandex.keys()])+"\n")
        self.plan=plandex
        self.html_title="<a href='"+self.title+".html' title='"+self.summary.replace("\n"," ")+"'>"+self.title+"</a>"
        self.content=txt

    def __str__(self):
        return self.content

def genelookup_factory(proteomics):
    ws = pd.read_excel(proteomics, index_col=0)
    nws=ws/ws.median()
    def genelookup_fun(gene):
        try:
            rowslice=nws.loc[gene]
            groupings=DefaultOrderedDict(list)
            for col in rowslice.index:
                #print("what you need to still do is group the .x together and make a json string for the googledex and sub the $google, also change the genedex")
                groupings[re.match("(.*?)\.",col).group(1)].append(rowslice.loc[col])
            m=max([len(groupings[x]) for x in groupings])
            for x in groupings:  # dict comprehension not a good idea as it is a DefaultOrderedDict
                groupings[x]+=[None for y in range(m-len(groupings[x]))]
            #return json.dumps([[g]+groupings[g] for g in groupings])
            return json.dumps([{"name":g, "type":'box',"y":groupings[g], "boxpoints": 'all'} for g in groupings])
        except KeyError:
            print("\t No data: ", gene)
            return None
    return genelookup_fun

def staticBioWikiMaker(markdown, proteomics,entrytemplate,categorytemplate, indexpage):
    #make the database
    mdtxt=open(markdown,encoding="utf-8").read()
    db=[mdentry(a) for a in re.split("\n(?=\# )",mdtxt)]  #fix again absent spaces and firt wonky one TODO
    modmould=Template(open(entrytemplate).read())
    #figure out all categories
    protocatdex=defaultdict(list)
    for topic in db:
        for tag in topic.tags:
            protocatdex[tag].append(topic.html_title)
    catmould=Template(open(categorytemplate).read())
    catdex={}
    for tag in protocatdex:
        catdex[tag]=catmould.substitute(title=tag, content=", ".join(protocatdex[tag]))
    menu=' '.join(['<li>'+topic.html_title+'</li>' for topic in db])
    genelookup=genelookup_factory(proteomics)
    # make html!
    for topic in db:
        html=md.markdown(topic.content, ['markdown.extensions.extra']).replace("<img ","<img width=100% ").replace("<table>","<table class='table table-striped' >") #fix
        html=modmould.safe_substitute(title=topic.title, content=html, menu=menu)
        html=Template(html).safe_substitute(catdex)
        #get the plan data.
        if topic.plan:
            plandex={strategy: catmould.substitute(title=strategy, content=", ".join(topic.plan[strategy])) for strategy in topic.plan}
            html=Template(html).safe_substitute(plandex)
        #get the gene data.
        if topic.genes:
            genedatadex={gene: genelookup(gene) for gene in topic.genes}
            jsdex={gene: "Plotly.newPlot('chart_"+gene+"', "+genedatadex[gene]+", {yaxis: {title: 'Spectral counts relative to median'}});" for gene in genedatadex if genedatadex[gene]}
            js="\n".join(jsdex.values())
            html=Template(html).safe_substitute(js=js)
            content={}
            chart={gene: '<div id="chart_'+gene+'"></div>' for gene in jsdex}
            for gene in topic.genes:
                content[gene]="<a href='http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&redirect=T&object="+gene+"'>Ecocyc</a><br />"
                if gene not in chart:
                    content[gene]+=""
                else:
                    content[gene]+=chart[gene]
            genedex={gene: catmould.substitute(title=gene, content=content[gene]) for gene in topic.genes}
            html=Template(html).safe_substitute(genedex)
        open("site/"+topic.title+".html","w", encoding="utf-8").write(html)
    ###special pages
    #main. for now I'll botch the conversion of SVG to links I am not well versed with ET. Please do not copy.
    #svg=open("data/images/map.svg").read()
    #svg=re.sub("<\!.*?>","",svg)
    #svgtree=et.fromstring(svg)
    et.register_namespace("","http://www.w3.org/2000/svg")
    svgtree=et.parse("data/map.svg").getroot()
    #<a xlink:href="#target">
    svgtree.attrib["width"]="100%"
    del svgtree.attrib["height"]
    for group in svgtree:
        if group.attrib["id"]=='clickareas':
            elball=[elem for elem in group] #Neeeded
            for elem in elball:
                if "id" in elem.attrib:
                    if 1==0:
                        group.remove(elem)
                        mod=elem.attrib["id"].replace("_x5F_","_")
                        for x in db: #a curious need for a if mod in [x.title for x in db] would lose the x.
                            if mod == x.title:
                                summary=x.summary
                                break
                        else:
                            summary=mod
                            print("Could not find ", mod, " among ", db)
                        n=et.SubElement(group,"a",{"xlink:href":mod+".html","xlink:title":summary})
                        n.insert(0,elem)
                    else:
                        elem.attrib["class"]="figlink"
                        mod=elem.attrib["id"].replace("_x5F_","_")
                        for x in db: #a curious need for a if mod in [x.title for x in db] would lose the x.
                            if mod == x.title:
                                summary=x.summary
                                break
                        else:
                            summary=mod
                            print("Could not find ", mod, " among ", db)
                        elem.attrib["xlink:href"]=mod+".html"
                        elem.attrib["xlink:title"]=summary
                        elem.attrib["onclick"]="javascript:location.href='"+mod+".html'"
    index=Template(open(indexpage).read()).substitute(title="Gene network", content=et.tostring(svgtree, encoding="utf-8").decode(), menu=menu)
    open("site/index.html","w", encoding="utf-8").write(index)
    #changes.
    changedex={x:set() for x in mdentry.plan_options}
    for topic in db:
        for strategy in topic.plan:
            changedex[strategy].update(topic.plan[strategy])
    changetxt="<div class='alert alert-danger' role='alert'>This section is deprecated.</div>\n"+"\n".join(["# "+strategy.capitalize()+"\n"+"\n".join(changedex[strategy]) for strategy in mdentry.plan_options])
    open("site/changes.html","w", encoding="utf-8").write(modmould.safe_substitute(menu=menu, title="changes",content=md.markdown(changetxt, ['markdown.extensions.extra'])))

if __name__=="__main__":
    staticBioWikiMaker("modules.md",'proteomics.xlsx',"data/template.html", "data/cat_template.html","data/itemplate.html")
