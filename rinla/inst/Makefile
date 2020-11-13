NEWS.pdf : NEWS.Rd
	R CMD Rd2pdf  -o NEWS.pdf NEWS.Rd

NEWS.md: NEWS.Rd Makefile
	cat NEWS.Rd \
	| sed "s/\\\newcommand{.*//" \
	| sed "s/\\\name{.*//" \
	| sed "s/\\\title{.*//" \
	| sed "s/\\\encoding{.*//" \
	| sed "s/.*\\\section{\([^}]*\)}{.*/\\# \1/" \
	| sed "s/.*\\\section{\([^}]*\)}{.*/\\# \1/" \
	| sed "s/.*\\\section{\([^}]*\)}{.*/\\# \1/" \
	| sed "s/CHANGES IN VERSION xx.xx.xx/INLA 99.99.99/" \
	| sed "s/CHANGES IN VERSION/INLA/" \
	| sed "s/.*\\\subsection{\([^}]*\)}{.*/\\## \1/" \
	| sed "s/.*\\\itemize{.*//" \
	| sed "s/.*\\\item/* /" \
	| sed "s/\\\code{\([^ }]*\)}/\`\1\`/g" \
	| sed "s/^[ ]*}//" \
	> NEWS.md
	@echo "NEWS.md generated from NEWS.rd."
	@echo "Check the contents. If OK, then move NEWS.md"
	@echo "to the package root and remove inst/NEWS.Rd (and inst/Makefile)"
	@echo 'The R news("INLA") output can first be checked with'
	@echo '  {'
	@echo '    db <- tools:::.build_news_db_from_package_NEWS_md("NEWS.md");'
	@echo '    attr(db, "package") <- "INLA";'
	@echo '    db'
	@echo '  }'
