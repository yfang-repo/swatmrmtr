# Import variables/options/rules from PETSc.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
# # Those files also include 
# # ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables
# # and 
# # ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscrules
# # so that PFLOTRAN will be built with the same options as 
# # the petsc configured in $PETSC_DIR/$PETSC_ARCH
#
#MYFLAGS = -I.
FCPPFLAGS += -I${chem_dir}
#FC_FLAGS = -fPIC -g -CB
FC_FLAGS = -fPIC -O2
#FLINKER = mpif90 -fPIC -g -CB -O2
#FLINKER = mpif90 -fPIC -CB -O2
FLINKER = mpif90 -fPIC -CB -g
chem_dir = ../pflotran/src/pflotran
LIBS += -L${chem_dir} -lpflotranchem

###############################################################################
# Assign additional compiler/preprocessor flags
# ###############################################################################
#
# # These flags are supplemental to the PETSc flags
CFLAGS   =
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}
#
#
swat_obj = zeroini.o \
	zero_urbn.o \
	zero2.o \
	zero1.o \
	zero0.o \
	ysed.o \
	xmon.o \
	writem.o \
	writed.o \
	writeaa.o \
	writea.o \
	wndgen.o \
	wmeas.o \
	wetlan.o \
	weatgn.o \
	watuse.o \
	wattable.o \
	water_hru.o \
	watbal.o \
	washp.o \
	volq.o \
	virtual.o \
	vbl.o \
	varinit.o \
	urbanhr.o \
	urban.o \
	urb_bmp.o \
	ttcoef_wway.o \
	ttcoef.o \
	tstr.o \
	transfer.o \
	tran.o \
	tmeas.o \
	tillmix.o \
	tillfactor.o \
	theta.o \
	tgen.o \
	tair.o \
	swu.o \
	sweep.o \
	swbl.o \
	surq_greenampt.o \
	surq_daycn.o \
	surfstor.o \
	surfst_h2o.o \
	surface.o \
	sumv.o \
	sumhyd.o \
	subyr.o \
	subwq.o \
	substor.o \
	submon.o \
	subday.o \
	subbasin.o \
	subaa.o \
	sub_subbasin.o \
	structure.o \
	storeinitial.o \
	stdaa.o \
	std3.o \
	std2.o \
	std1.o \
	solt.o \
	solp.o \
	soil_write.o \
	soil_phys.o \
	soil_par.o \
	soil_chem.o \
	snom.o \
	smeas.o \
	slrgen.o \
	simulate.o \
	sim_inityr.o \
	sim_initday.o \
	schedule_ops.o \
	sched_mgt.o \
	saveconc.o \
	save.o \
	sat_excess.o \
	rtsed_yangsand.o \
	rtsed_kodatie.o \
	rtsed_bagnold.o \
	rtsed_Molinas_Wu.o \
	rtsed.o \
	rtpest.o \
	rtout.o \
	rtmusk.o \
	rthsed.o \
	rthr.o \
	rthpest.o \
	rthmusk.o \
	rteinit.o \
	rtday.o \
	rtbact.o \
	rsedyr.o \
	rsedmon.o \
	rseday.o \
	rsedaa.o \
	routres.o \
	routeunit.o \
	routels.o \
	route.o \
	rootfr.o \
	rhgen.o \
	rewind_init.o \
	resnut.o \
	resinit.o \
	reshr.o \
	resetlu.o \
	resbact.o \
	res.o \
	regres.o \
	recyear.o \
	recmon.o \
	rechour.o \
	recday.o \
	reccnst.o \
	readyr.o \
	readwwq.o \
	readwus.o \
	readwgn.o \
	readurban.o \
	readtill.o \
	readswq.o \
	readsub.o \
	readsol.o \
	readsno.o \
	readseptwq.o \
	readsepticbz.o \
	readsdr.o \
	readru.o \
	readrte.o \
	readres.o \
	readpnd.o \
	readplant.o \
	readpest.o \
	readops.o \
	readmon.o \
	readmgt.o \
	readlwq.o \
	readlup.o \
	readinpt.o \
	readhru.o \
	readgw.o \
	readfile.o \
	readfig.o \
	readfert.o \
	readfcst.o \
	readcnst.o \
	readchm.o \
	readbsn.o \
	readatmodep.o \
	reachout.o \
	rchyr.o \
	rchuse.o \
	rchmon.o \
	rchinit.o \
	rchday.o \
	rchaa.o \
	ran1.o \
	qman.o \
	psed.o \
	print_hyd.o \
	potholehr.o \
	pothole.o \
	pondhr.o \
	pond.o \
	pminrl2.o \
	pminrl.o \
	pmeas.o \
	plantop.o \
	plantmod.o \
	pkq.o \
	pgenhr.o \
	pgen.o \
	pesty.o \
	pestw.o \
	pestlch.o \
	percmicro.o \
	percmain.o \
	percmacro.o \
	ovr_sed.o \
	origtile.o \
	orgncswat.o \
	orgn.o \
	operatn.o \
	openwth.o \
	nuts.o \
	nup.o \
	nrain.o \
	npup.o \
	noqual.o \
	nminrl.o \
	nlch.o \
	nitvol.o \
	nfix.o \
	newtillmix.o \
	ndenit.o \
	modparm.o \
	lwqdef.o \
	lids.o \
	lidinit.o \
	lid_raingarden.o \
	lid_porpavement.o \
	lid_greenroof.o \
	lid_cistern.o \
	layersplit.o \
	latsed.o \
	lakeq.o \
	killop.o \
	jdt.o \
	irrsub.o \
	irrigate.o \
	irr_res.o \
	irr_rch.o \
	impndyr.o \
	impndmon.o \
	impndday.o \
	impndaa.o \
	impnd_init.o \
	icl.o \
	hydroinit.o \
	hruyr.o \
	hrupondhr.o \
	hrupond.o \
	hrumon.o \
	hruallo.o \
	hruaa.o \
	hmeas.o \
	hhnoqual.o \
	headout.o \
	header.o \
	harvkillop.o \
	harvgrainop.o \
	harvestop.o \
	h2omgt_init.o \
	gwnutr.o \
	gwmod_deep.o \
	gwmod.o \
	gw_no3.o \
	grow.o \
	graze.o \
	grass_wway.o \
	getallo.o \
	gcycl.o \
	finalbal.o \
	filtw.o \
	filter.o \
	fert.o \
	expo.o \
	etpot.o \
	etact.o \
	estimate_ksat.o \
	erfc.o \
	enrsb.o \
	eiusle.o \
	ee.o \
	dstn1.o \
	drains.o \
	dormant.o \
	distrib_bmps.o \
	depstor.o \
	decay.o \
	dailycn.o \
	curno.o \
	crackvol.o \
	crackflow.o \
	confert.o \
	conapply.o \
	command.o \
	clicon.o \
	clgen.o \
	chkcst.o \
	cfactor.o \
	carbon_new.o \
	caps.o \
	canopyint.o \
	burnop.o \
	buffer.o \
	bmpinit.o \
	bmpfixed.o \
	bmp_wet_pond.o \
	bmp_sed_pond.o \
	bmp_sand_filter.o \
	bmp_ri_pond.o \
	bmp_det_pond.o \
	biozone.o \
	biofilm.o \
	bacteria.o \
	aveval.o \
	autoirr.o \
	aunif.o \
	atri.o \
	ascrv.o \
	apply.o \
	apex_day.o \
	anfert.o \
	alph.o \
	allocate_parms.o \
	albedo.o \
	addh.o \
	hruday.o \
	carbon_zhang2.o \
	NCsed_leach.o


main.o: main.F90
	${FLINKER} -c -o $@ $^ ${FCPPFLAGS} ${PETSC_LIB} ${LIBS}

hhwatqual.o: hhwatqual.f
	${FLINKER} -c -o $@ $^ ${FCPPFLAGS} ${PETSC_LIB} ${LIBS}
watqual.o: watqual.f
	${FLINKER} -c -o $@ $^ ${FCPPFLAGS} ${PETSC_LIB} ${LIBS}
watqual2.o: watqual2.f
	${FLINKER} -c -o $@ $^ ${FCPPFLAGS} ${PETSC_LIB} ${LIBS}

pflotran_rxn_obj = pflotran_rxn.o
pflotran_rxn.o: pflotran_rxn.F90
	${FLINKER} -c -o $@ $^ ${FCPPFLAGS} ${PETSC_LIB} ${LIBS}

swatmrmtr.e : $(pflotran_rxn_obj) hhwatqual.o watqual2.o watqual.o main.o $(swat_obj)
	${FLINKER} -o swatmrmtr.e ${swat_obj} hhwatqual.o watqual.o watqual2.o main.o ${pflotran_rxn_obj} ${PETSC_LIB} ${LIBS}

###############################################################################
# Dependencies
# # Dependencies stemming from "use" statements.
# # These ensure that the module files are built in the correct order.
# # The source file (i.e. file.F90) is intentionally omitted; 
# # it is picked up by the default/implicit build rule (defined by PETSc)
# ###############################################################################
#
addh.o:      modparm.o
albedo.o:      modparm.o
allocate_parms.o:      modparm.o
alph.o:      modparm.o
anfert.o:      modparm.o
apex_day.o:      modparm.o
apply.o:      modparm.o
autoirr.o:      modparm.o
bacteria.o:      modparm.o
biofilm.o:      modparm.o
biozone.o:      modparm.o
bmp_det_pond.o:      modparm.o
bmpfixed.o:      modparm.o
bmpinit.o:      modparm.o
bmp_ri_pond.o:      modparm.o
bmp_sand_filter.o:      modparm.o
bmp_sed_pond.o:      modparm.o
bmp_wet_pond.o:      modparm.o
bmp_wet_pond.o:      modparm.o
buffer.o:      modparm.o
burnop.o:      modparm.o
canopyint.o:      modparm.o
carbon_new.o:      modparm.o
cfactor.o:      modparm.o
clgen.o:      modparm.o
clicon.o:      modparm.o
command.o:      modparm.o
conapply.o:      modparm.o
confert.o:      modparm.o
crackflow.o:      modparm.o
crackvol.o:      modparm.o
curno.o:      modparm.o
dailycn.o:      modparm.o
decay.o:      modparm.o
depstor.o:      modparm.o
distrib_bmps.o:      modparm.o
dormant.o:      modparm.o
drains.o:      modparm.o
eiusle.o:      modparm.o
enrsb.o:      modparm.o
etact.o:      modparm.o
etpot.o:      modparm.o
fert.o:      modparm.o
filter.o:      modparm.o
filtw.o:      modparm.o
finalbal.o:      modparm.o
gcycl.o:      modparm.o
getallo.o:      modparm.o
grass_wway.o:      modparm.o
graze.o:      modparm.o
grow.o:      modparm.o
gwmod_deep.o:      modparm.o
gwmod.o:      modparm.o
gw_no3.o:      modparm.o
gwnutr.o:      modparm.o
h2omgt_init.o:      modparm.o
harvestop.o:      modparm.o
harvgrainop.o:      modparm.o
harvkillop.o:      modparm.o
header.o:      modparm.o
headout.o:      modparm.o
hhnoqual.o:      modparm.o
hhwatqual.o:      modparm.o
hmeas.o:      modparm.o
hruaa.o:      modparm.o
hruallo.o:      modparm.o
hrumon.o:      modparm.o
hrupond.o:      modparm.o
hrupondhr.o:      modparm.o
hruyr.o:      modparm.o
hydroinit.o:      modparm.o
icl.o:      modparm.o
impndaa.o:      modparm.o
impndday.o:      modparm.o
impnd_init.o:      modparm.o
impndmon.o:      modparm.o
impndyr.o:      modparm.o
irrigate.o:      modparm.o
irr_rch.o:      modparm.o
irr_res.o:      modparm.o
irrsub.o:      modparm.o
killop.o:      modparm.o
lakeq.o:      modparm.o
latsed.o:      modparm.o
layersplit.o:      modparm.o
lid_cistern.o:      modparm.o
lid_greenroof.o:      modparm.o
lidinit.o:      modparm.o
lid_porpavement.o:      modparm.o
lid_raingarden.o:      modparm.o
lids.o:      modparm.o
lwqdef.o:      modparm.o
main.o:      modparm.o $(chem_dir)/libpflotranchem.a
ndenit.o:	modparm.o
newtillmix.o:      modparm.o
nfix.o:      modparm.o
nitvol.o:      modparm.o
nlch.o:      modparm.o
nminrl.o:      modparm.o
noqual.o:      modparm.o
npup.o:      modparm.o
nrain.o:      modparm.o
nup.o:      modparm.o
openwth.o:      modparm.o
operatn.o:      modparm.o
orgncswat.o:      modparm.o
orgn.o:      modparm.o
origtile.o:      modparm.o
ovr_sed.o:      modparm.o
percmacro.o:      modparm.o
percmain.o:      modparm.o
percmicro.o:      modparm.o
pestlch.o:      modparm.o
pestw.o:      modparm.o
pesty.o:      modparm.o
pgen.o:      modparm.o
pgenhr.o:      modparm.o
pkq.o:      modparm.o
plantmod.o:      modparm.o
plantop.o:      modparm.o
pmeas.o:      modparm.o
pminrl2.o:      modparm.o
pminrl.o:      modparm.o
pond.o:      modparm.o
pondhr.o:      modparm.o
pothole.o:      modparm.o
potholehr.o:	  modparm.o
print_hyd.o:      modparm.o
psed.o:      modparm.o
rchaa.o:      modparm.o
rchday.o:      modparm.o
rchinit.o:      modparm.o
rchmon.o:      modparm.o
rchuse.o:      modparm.o
rchyr.o:      modparm.o
reachout.o:      modparm.o
readatmodep.o:      modparm.o
readbsn.o:      modparm.o
readchm.o:      modparm.o
readcnst.o:      modparm.o
readfcst.o:      modparm.o
readfert.o:      modparm.o
readfig.o:      modparm.o
readfile.o:      modparm.o
readgw.o:      modparm.o
readhru.o:      modparm.o
readinpt.o:      modparm.o
readlup.o:      modparm.o
readlwq.o:      modparm.o
readmgt.o:      modparm.o
readmon.o:      modparm.o
readops.o:      modparm.o
readpest.o:      modparm.o
readplant.o:      modparm.o
readpnd.o:      modparm.o
readres.o:      modparm.o
readrte.o:      modparm.o
readru.o:      modparm.o
readsdr.o:      modparm.o
readsepticbz.o:      modparm.o
readseptwq.o:      modparm.o
readsno.o:      modparm.o
readsol.o:      modparm.o
readsub.o:      modparm.o
readswq.o:      modparm.o
readtill.o:      modparm.o
readurban.o:      modparm.o
readwgn.o:      modparm.o
readwus.o:      modparm.o
readwwq.o:      modparm.o
readyr.o:      modparm.o
reccnst.o:      modparm.o
recday.o:      modparm.o
rechour.o:      modparm.o
recmon.o:      modparm.o
recyear.o:      modparm.o
regres.o:      modparm.o
resbact.o:      modparm.o
resetlu.o:      modparm.o
res.o:      modparm.o
reshr.o:      modparm.o
resinit.o:      modparm.o
resnut.o:      modparm.o
rewind_init.o:      modparm.o
rhgen.o:      modparm.o
rootfr.o:	modparm.o
route.o:      modparm.o
routels.o:      modparm.o
routeunit.o:      modparm.o
routres.o:      modparm.o
rsedaa.o:      modparm.o
rseday.o:      modparm.o
rsedmon.o:      modparm.o
rsedyr.o:      modparm.o
rtbact.o:      modparm.o
rtday.o:      modparm.o
rteinit.o:      modparm.o
rthmusk.o:      modparm.o
rthpest.o:      modparm.o
rthr.o:      modparm.o
rthsed.o:	modparm.o
rtmusk.o:      modparm.o
rtout.o:      modparm.o
rtpest.o:      modparm.o
rtsed_bagnold.o:      modparm.o
rtsed.o:      modparm.o
rtsed_kodatie.o:      modparm.o
rtsed_Molinas_Wu.o:      modparm.o
rtsed_yangsand.o:      modparm.o
sat_excess.o:      modparm.o
saveconc.o:      modparm.o
save.o:      modparm.o
sched_mgt.o:      modparm.o
schedule_ops.o:      modparm.o
sim_initday.o:      modparm.o
sim_inityr.o:      modparm.o
simulate.o:      modparm.o
slrgen.o:      modparm.o
smeas.o:      modparm.o
snom.o:      modparm.o
soil_chem.o:      modparm.o
soil_par.o:      modparm.o
soil_phys.o:      modparm.o
soil_write.o:      modparm.o
solp.o:      modparm.o
solt.o:      modparm.o
std1.o:      modparm.o
std2.o:      modparm.o
std3.o:      modparm.o
stdaa.o:      modparm.o
storeinitial.o:      modparm.o
structure.o:      modparm.o
subaa.o:      modparm.o
subbasin.o:      modparm.o
subday.o:      modparm.o
submon.o:      modparm.o
substor.o:      modparm.o
sub_subbasin.o:      modparm.o
subwq.o:      modparm.o
subyr.o:      modparm.o
sumhyd.o:      modparm.o
sumv.o:      modparm.o
surface.o:      modparm.o
surfst_h2o.o:      modparm.o
surfstor.o:      modparm.o
surq_daycn.o:      modparm.o
surq_greenampt.o:      modparm.o
swbl.o:      modparm.o
sweep.o:      modparm.o
swu.o:      modparm.o
tair.o:      modparm.o
tgen.o:      modparm.o
tillfactor.o:      modparm.o
tillmix.o:      modparm.o
tmeas.o:      modparm.o
tran.o:      modparm.o
transfer.o:      modparm.o
tstr.o:      modparm.o
ttcoef.o:      modparm.o
ttcoef_wway.o:      modparm.o
urban.o:      modparm.o
urbanhr.o:      modparm.o
urb_bmp.o:      modparm.o
varinit.o:      modparm.o
vbl.o:      modparm.o
virtual.o:      modparm.o
volq.o:      modparm.o
washp.o:      modparm.o
watbal.o:      modparm.o
water_hru.o:      modparm.o
watqual2.o:      modparm.o pflotran_rxn.o
watqual.o:      modparm.o pflotran_rxn.o
wattable.o:      modparm.o
watuse.o:      modparm.o
weatgn.o:      modparm.o
wetlan.o:      modparm.o
wmeas.o:      modparm.o
wndgen.o:      modparm.o
writeaa.o:      modparm.o
writea.o:      modparm.o
writed.o:      modparm.o
writem.o:      modparm.o
xmon.o:      modparm.o
ysed.o:      modparm.o
zero0.o:      modparm.o
zero1.o:      modparm.o
zero2.o:      modparm.o
zeroini.o:      modparm.o
zero_urbn.o:      modparm.o
carbon_zhang2.o:        modparm.o
hruday.o:      modparm.o
NCsed_leach.o:      modparm.o
pflotran_rxn.o: $(chem_dir)/libpflotranchem.a modparm.o


echo-flags :
	@echo ; \
	echo "PFLOTRAN source"; \
	echo "SRC_DIR = $(SRC_DIR)"; \
	echo ; \
	echo "PETSc environment"; \
	echo "PETSC_DIR = $(PETSC_DIR)"; \
	echo "PETSC_ARCH = $(PETSC_ARCH)"; \
	echo ; \
	echo "Flags defined in this makefile"; \
	echo "CFLAGS = $(CFLAGS)"; \
	echo "FFLAGS = $(FFLAGS)"; \
	echo ; \
	echo "MYFLAGS = $(MYFLAGS)"; \
	echo "CPPFLAGS = $(CPPFLAGS)"; \
	echo "FPPFLAGS = $(FPPFLAGS)"; \
	echo ; \
	echo "Flags invoked when compiling fortran source"; \
	echo "FC_FLAGS = $(FC_FLAGS)"; \
	echo "FFLAGS = $(FFLAGS)"; \
	echo "FCPPFLAGS =  $(FCPPFLAGS)"; \
	echo ; \
	echo "Included Libraries" ;\
	echo "LIBS = $(LIBS)"; \
	echo ; \
	echo "PETSC_LIB = $(PETSC_LIB)"; \
	echo ; \
	echo "PETSC_LIB = $(FLINKER)"; \
	echo ; \
	echo "chem_LIB = $(chem_lib_path)"; \
	echo ; \

