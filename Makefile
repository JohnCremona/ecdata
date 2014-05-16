#Makefile for ECDATA directory

CC = gcc

FTP_HOST = warwick
FTP_DIR = public_html/ftp/data
# FTP_DIR = ./ftp/data

ALLCURVES = allcurves.?0000-?9999 allcurves.??0000-??9999
APLIST = aplist.?0000-?9999 aplist.??0000-??9999
BIGSHA = allbigsha.?0000-?9999 allbigsha.??0000-??9999
COUNT = count.?0000-?9999 count.??0000-??9999
CURVES = curves.?0000-?9999 curves.??0000-??9999
DEGPHI = degphi.?0000-?9999 degphi.??0000-??9999
ALLDEGPHI = alldegphi.?0000-?9999 alldegphi.??0000-??9999
GENS =
ALLGENS = allgens.?0000-?9999 allgens.??0000-??9999
ALLISOG = allisog.?0000-?9999 allisog.??0000-??9999
BSD = bsd.1-1000 allbsd.?0000-?9999 allbsd.??0000-??9999
PARICURVES = paricurves.?0000-?9999 paricurves.??0000-??9999
INTPTS = intpts.?0000-?9999 intpts.??0000-??9999
GALREPS = galrep.?0000-?9999 galrep.??0000-??9999
OPTIMAL = optimality.??

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
optimal: $(OPTIMAL)
	@echo $(OPTIMAL)

HTMLDATAFILES = shas.html table.html curves.1-1000.html
TEXTFILES = manin.txt INDEX.html release_notes.txt
DATAFILES =  $(ALLCURVES) $(APLIST) $(BIGSHA) $(COUNT) $(DEGPHI) $(ALLDEGPHI) $(ALLGENS) $(BSD) $(ALLISOG) $(PARICURVES) $(INTPTS) $(GALREPS) $(OPTIMAL)
FTPFILES = $(DATAFILES) $(TEXTFILES) $(HTMLDATAFILES)

ftp:  $(FTPFILES)
	for f in $(DATAFILES); \
	do \
	     gzip -c $${f} > ftpdir/$${f}.gz; \
	done;
	for f in $(TEXTFILES) $(HTMLDATAFILES); \
	do \
	     rsync -av $${f} ftpdir/$${f}; \
	done; \
        echo Updating $(FTP_HOST):$(FTP_DIR); \
        rsync -avz ftpdir/ $(FTP_HOST):$(FTP_DIR)/; \
	ssh warwick chmod a+rx $(FTP_DIR)
	ssh warwick chmod a+r $(FTP_DIR)/*

tar_old: $(FTPFILES)
	rm -f ftpdata* ftpfiles
	touch ftpfiles
	for f in $(FTPFILES); \
	do echo $${f} >> ftpfiles; done
	tar -zcf ftpdata.tgz --files-from=ftpfiles
	mv ftpdata.tgz $(FTP_DIR)/..
	chmod 644 $(FTP_DIR)/../ftpdata.tgz

DATE = $(shell date +%Y-%m-%d )
tar: $(FTPFILES)
	rm -f ecdata/*
	for f in $(FTPFILES); do ln -s $(PWD)/$${f} ecdata/$${f}; done
	tar -zchf ecdata-$(DATE).tgz ecdata
	scp ecdata-$(DATE).tgz $(FTP_HOST):$(FTP_DIR)/..
	ssh $(FTP_HOST) chmod a+r $(FTP_DIR)/../ecdata-$(DATE).tgz
