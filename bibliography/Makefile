all: README.md

biblio.png: biblio.bib
	@echo "\033[1m>\033[0m \033[34m$@\033[0m"
	./authorsgraph.py $^ | neato -T png -o $@

README.md: biblio.bib pdf biblio.png
	@echo "\033[1m>\033[0m \033[34m$@\033[0m"
	./bibdisplay.py > $@

clean:
	@echo "\033[1m>\033[0m \033[31m$@\033[0m"
	rm README.md biblio.png
