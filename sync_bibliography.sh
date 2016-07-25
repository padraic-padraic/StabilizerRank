#!/usr/bin/env bash
bibfile="/Users/padraic/Documents/Bibtex/MResProject.bib"
mybib="/Users/padraic/Documents/UCL/MRes/Project/MResProject.bib"
if [  -f ${bibfile} ]; then
    cp -f $bibfile $mybib
fi

