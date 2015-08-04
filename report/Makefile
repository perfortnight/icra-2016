all: paper

paper: root.tex references.bib
	pdflatex -jobname=icra2016 -file-line-error root.tex
	bibtex icra2016.aux
	pdflatex -jobname=icra2016 -file-line-error root.tex
	pdflatex -jobname=icra2016 -file-line-error root.tex
	evince icra2016.pdf


clean:
	rm icra2016.aux icra2016.bbl icra2016.blg icra2016.log icra2016.out
