all: doc

doc:
	pdflatex main.tex ; pdflatex main.tex

clean:
	find . -type f \( -name "*~" -or -name "#*#" \
-or -name "*.pyc" -or -name "*.out" -or -name "*.aux" -or -name "*.bbl" \
-or -name "*.log" -or -name "*.toc" -or -name "*.blg" -or -name "*.dvi" \
-or -name "*-blx.bib" -or -name "*.run.xml" \) -delete
