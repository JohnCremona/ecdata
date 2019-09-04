#Makefile for ECDATA directory

CC = gcc

FTP_HOST = warwick
FTP_DIR = public_html/ftp/data
DIST_DIR = /home/jec/ecdata/ecdata

ALLCURVES = allcurves/allcurves.?0000-?9999 allcurves/allcurves.??0000-??9999
APLIST = aplist/aplist.?0000-?9999 aplist/aplist.??0000-??9999
BIGSHA = allbigsha/allbigsha.?0000-?9999 allbigsha/allbigsha.??0000-??9999
COUNT = count/count.?0000-?9999 count/count.??0000-??9999
CURVES = curves/curves.?0000-?9999 curves/curves.??0000-??9999
DEGPHI = degphi/degphi.?0000-?9999 degphi/degphi.??0000-??9999
ALLDEGPHI = alldegphi/alldegphi.?0000-?9999 alldegphi/alldegphi.??0000-??9999
GENS =
ALLGENS = allgens/allgens.?0000-?9999 allgens/allgens.??0000-??9999
ALLISOG = allisog/allisog.?0000-?9999 allisog/allisog.??0000-??9999
BSD = allbsd/bsd.1-1000 allbsd/allbsd.?0000-?9999 allbsd/allbsd.??0000-??9999
PARICURVES = paricurves/paricurves.?0000-?9999 paricurves/paricurves.??0000-??9999
INTPTS = intpts/intpts.?0000-?9999 intpts/intpts.??0000-??9999
GALREPS = galrep/galrep.?0000-?9999 galrep/galrep.??0000-??9999
TWOADIC = 2adic/2adic.?0000-?9999 2adic/2adic.??0000-??9999
OPTIMAL = opt_man/opt_man.?0000-?9999 opt_man/opt_man.??0000-??9999
GROWTH = growth/growth?d.*000-*999

allcurves: $(ALLCURVES)
	@echo $(ALLCURVES)
aplist: $(APLIST)
	@echo $(APLIST)
allbigsha: $(ALLBIGSHA)
	@echo $(ALLBIGSHA)
count: $(COUNT)
	@echo $(COUNT)
curves: $(CURVES)
	@echo $(CURVES)
degphi: $(DEGPHI)
	@echo $(DEGPHI)
alldegphi: $(ALLDEGPHI)
	@echo $(ALLDEGPHI)
allgens: $(ALLGENS)
	@echo $(ALLGENS)
allisog: $(ALLISOG)
	@echo $(ALLISOG)
bsd: $(BSD)
	@echo $(BSD)
paricurves: $(PARICURVES)
	@echo $(PARICURVES)
intpts: $(INTPTS)
	@echo $(INTPTS)
galreps: $(GALREPS)
	@echo $(GALREPS)
twoadic: $(TWOADIC)
	@echo $(TWOADIC)
optimal: $(OPTIMAL)
	@echo $(OPTIMAL)
growth: $(GROWTH)
	@echo $(GROWTH)

HTMLFILES = html/index.html html/shas.html html/table.html html/curves.1-1000.html
TEXTFILES = doc/manin.txt doc/file-format.txt doc/release_notes.md doc/merging.txt
DATAFILES =  $(ALLCURVES) $(APLIST) $(BIGSHA) $(COUNT) $(DEGPHI) $(ALLDEGPHI) $(ALLGENS) $(BSD) $(ALLISOG) $(PARICURVES) $(INTPTS) $(GALREPS) $(TWOADIC) $(OPTIMAL) $(GROWTH)
FTPFILES = $(DATAFILES) $(TEXTFILES) $(HTMLFILES)
DATASUBDIRS = allcurves aplist allbigsha count curves degphi alldegphi allgens allisog allbsd paricurves intpts galrep 2adic growth opt_man

commit: $(FTPFILES)
	git add $(DATAFILES)
	-git commit -m "updated data files"
	git add $(TEXTFILES)
	-git commit -m "updated text files (in master:./doc)"
	git add $(HTMLFILES)
	-git commit -m "updated html files (in master:./html)"
	git push
	git checkout gh-pages
	git checkout master:./html/ ./
	git add *.html
	-git commit -m "updated html files in gh-pages branch"
	git push
	git checkout master


DATE = $(shell date +%Y-%m-%d )
tar: $(FTPFILES)
	rm -rf $(DIST_DIR)
	mkdir -p $(DIST_DIR)
	for d in $(DATASUBDIRS) html doc; do mkdir -p $(DIST_DIR)/$${d}; done
	for f in $(FTPFILES); do ln -s $(PWD)/$${f} $(DIST_DIR)/$${f}; done
	cd $(DIST_DIR)/..
	tar -zchf ecdata-$(DATE).tgz ecdata
	scp ecdata-$(DATE).tgz $(FTP_HOST):$(FTP_DIR)/..
	ssh $(FTP_HOST) chmod a+r $(FTP_DIR)/../ecdata-$(DATE).tgz

# NB The Sage import script requires that the tarball extracts to a
# directory called 'ecdata'
