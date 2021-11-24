.PHONY: all viewpdf pdf clean

TARGET       = main
SOURCE_FILES = $(TARGET).tex $(wildcard */*.tex) $(wildcard */sections/*.tex)
BIB_FILES    = $(wildcard biblio/*.bib)
FIGURES      = $(wildcard */figures/*)

# Set the pdf reader according to the operating system
OS = $(shell uname)
ifeq ($(OS), Darwin)
	PDF_READER = open
endif
ifeq ($(OS), Linux)
	PDF_READER = evince
endif

all: pdf

viewpdf: pdf
	$(PDF_READER) $(TARGET).pdf &

pdf: $(TARGET).pdf

min: $(TARGET).min.pdf

%.min.pdf: %.pdf
	 gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dNOPAUSE -dQUIET -dBATCH -dPrinted=false -sOutputFile=$@ $<

one : $(SOURCE_FILES) $(BIB_FILES) $(FIGURES) these-dbl.cls
	pdflatex -jobname=$(TARGET) $(SOURCE_FILES)

$(TARGET).pdf: $(SOURCE_FILES) $(BIB_FILES) $(FIGURES) these-dbl.cls
	pdflatex $<
	biber $(TARGET)
	pdflatex $<
	pdflatex $<

clean:
	rm -f $(TARGET).{ps,pdf,min.pdf,bcf,run.xml}
	for suffix in dvi aux bbl blg toc ind out brf ilg idx synctex.gz log; do \
		find . -type d -name ".git" -prune -o -type f -name "*.$${suffix}" -print -exec rm {} \;  ; \
	done

