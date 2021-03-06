# Demo
Here is an [example version](https://rawgit.com/matteoferla/staticBioWikiMaker/master/site/index.html)

# staticBioWikiMaker
`staticBioWikiMaker(markdown, map, proteomics,entrytemplate,categorytemplate, indexpage)` makes a static set of webpages based on a markdown file.

The level 1 headers on the markdown file are splint into entries and inserted into a template (entrytemplate), which has a $title, $menu and $content keywords for substitution.

The level 2 headers are parse differently. Members should contain gene symbols prefixed with a dollar sign. The script will extract the proteomics data from the proteomics file and median normalise it and show on that entry the genes tagged in that manner.
Tags contains keywords with a dollar sigil that are common to other entries. The site will generate a navigation box for each keyword with links to all entries which mention it.

There are a few extra features not discussed.

Also to get google drive to serve the html as such the <id> needs to be put in https://googledrive.com/host/<id>

# Map

One major trick is to integrate the map with the pages.
The index page is an SVG with a layer called "clickareas" in front of the main picture and containing polygons of opacity 0% with ids matching the entries.
![Alt text](/demo/making_of_svg_IDs.jpg?raw=true "Making the interlink trick in Illustrator")

The script will link them and change the opacity onhover.
The entries will generate the dropdown menu `$menu`, but entries can be generated that are hardcoded in the template, such as an `about` page.
