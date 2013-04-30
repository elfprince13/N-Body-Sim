all: summary_abstract.bbl
	pdflatex summary_abstract.tex 
summary_abstract.bbl: summary_abstract.aux
	bibtex summary_abstract.aux

summary_abstract.aux:
	pdflatex summary_abstract.tex

clean:
	rm summary_abstract.pdf
	rm summary_abstract.log
	rm summary_abstract.bbl
	rm summary_abstract.blg
	rm summary_abstract.aux
