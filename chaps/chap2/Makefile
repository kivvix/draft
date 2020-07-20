MAIN_FILE = chap2
sections_files = $(wildcard sections/*.tex)

all: $(MAIN_FILE).pdf

min: $(MAIN_FILE).min.pdf

%.pdf: %.tex $(sections_files) biblio.bib
	echo $^
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

%.min.pdf: %.pdf
	 gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dNOPAUSE -dQUIET -dBATCH -dPrinted=false -sOutputFile=$@ $<

.PHONY: all clean clean_texshop mrproper

bib:
	biber -u -U -q chap2

clean: clean_texshop
	latexmk -c

clean_texshop:
	-rm *.run.xml

clean_biber:
	-rm -rf `biber --cache`

clean_min_pdf:
	-rm *.min.pdf

mrproper: clean_texshop clean_biber clean_min_pdf
	latexmk -C
