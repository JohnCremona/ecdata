#Makefile for ECDATA directory

CC = gcc

FTP_HOST = warwick
FTP_DIR = public_html/ftp/data
# FTP_DIR = ./ftp/data

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
OPTIMAL = optimality/optimality.??

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

HTMLDATAFILES = shas.html table.html curves.1-1000.html
TEXTFILES = manin.txt INDEX.html release_notes.txt
DATAFILES =  $(ALLCURVES) $(APLIST) $(BIGSHA) $(COUNT) $(DEGPHI) $(ALLDEGPHI) $(ALLGENS) $(BSD) $(ALLISOG) $(PARICURVES) $(INTPTS) $(GALREPS) $(TWOADIC) $(OPTIMAL)
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
	ssh warwick chmod -R a+rX $(FTP_DIR)

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
	rm -f ecdata/*.txt
	rm -f ecdata/*.html
	rm -f ecdata/*/*
	for f in $(FTPFILES); do ln -s $(PWD)/$${f} ecdata/$${f}; done
	tar -zchf ecdata-$(DATE).tgz ecdata
	scp ecdata-$(DATE).tgz $(FTP_HOST):$(FTP_DIR)/..
	ssh $(FTP_HOST) chmod a+r $(FTP_DIR)/../ecdata-$(DATE).tgz
