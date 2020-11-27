#Makefile for ECDATA directory

CC = gcc

FTP_DIR = $(HOME)/JohnCremona.github.io/ftp/data
DIST_DIR = $(HOME)/ecdata/ecdata

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
GROWTH = growth/*/growth*.*000-*999
IWASAWA = iwasawa/iwasawa.?0000-?9999 iwasawa/iwasawa.1?0000-1?9999

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
iwasawa: $(IWASAWA)
	@echo $(IWASAWA)

HTMLFILES = docs/index.html docs/shas.html docs/table.html docs/curves.1-1000.html
TEXTFILES = docs/manin.txt docs/file-format.txt docs/release_notes.md docs/merging.txt
DATAFILES =  $(ALLCURVES) $(APLIST) $(BIGSHA) $(COUNT) $(DEGPHI) $(ALLDEGPHI) $(ALLGENS) $(BSD) $(ALLISOG) $(PARICURVES) $(INTPTS) $(GALREPS) $(TWOADIC) $(OPTIMAL) $(GROWTH) $(IWASAWA)
FTPFILES = $(DATAFILES) $(TEXTFILES) $(HTMLFILES)
DATASUBDIRS = allcurves aplist allbigsha count curves degphi alldegphi allgens allisog allbsd paricurves intpts galrep 2adic growth opt_man iwasawa

commit: $(FTPFILES)
	git add $(DATAFILES)
	-git commit -m "updated data files"
	git add $(TEXTFILES)
	-git commit -m "updated text files (in master:./docs)"
	git add $(HTMLFILES)
	-git commit -m "updated html files (in master:./docs)"
	git push

DATE = $(shell date +%Y-%m-%d )
tar: $(FTPFILES)
	rm -rf $(DIST_DIR)
	mkdir -p $(DIST_DIR)
	for d in $(DATASUBDIRS) docs; do mkdir -p $(DIST_DIR)/$${d}; done
	for f in $(FTPFILES); do ln -s $(PWD)/$${f} $(DIST_DIR)/$${f}; done
	cd $(DIST_DIR)/..
	tar -zchf $(FTP_DIR)/ecdata-$(DATE).tgz ecdata
	rm -rf $(DIST_DIR)

# NB The Sage import script requires that the tarball extracts to a
# directory called 'ecdata'
