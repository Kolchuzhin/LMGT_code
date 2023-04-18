*version 6.0 3913105963
lib C:\MSIMEV60\LIB\abm.slb 0
lib C:\MSIMEV60\LIB\analog.slb 2746273205
lib C:\MSIMEV60\LIB\breakout.slb 3690318909
lib C:\MSIMEV60\LIB\eval.slb 0
lib C:\MSIMEV60\LIB\port.slb 0
lib C:\MSIMEV60\LIB\source.slb 2970098879
lib C:\MSIMEV60\LIB\special.slb 25865778
lib C:\MSIMEV60\LIB\7400.slb 0
lib C:\MSIMEV60\LIB\mechanik.slb 2050559946
lib C:\MSIMEV60\LIB\misc.slb 248682172
u 132
ABM1? 7
PM? 2
INTEG? 3
DIFF? 3
V? 5
R? 2
? 4
ABM? 2
@libraries
@analysis
.TRAN 1 0 0 0
+0 10ms
+1 10s
+3 10ms
.STEP 0 0 0
+ 0 V1
+ 4 0
+ 5 20000
+ 6 10
.PROBE 0 0 1 1 0 2
@targets
@attributes
@translators
a 0 u 13 0 0 0 hln 100 PSPICE=PSPICE
a 0 u 13 0 0 0 hln 100 PADS=PCB
a 0 u 13 0 0 0 hln 100 PCAD=PCB
a 0 u 13 0 0 0 hln 100 PROTEL=PCB
a 0 u 13 0 0 0 hln 100 TANGO=PCB
a 0 u 13 0 0 0 hln 100 CADSTAR=PCB
a 0 u 13 0 0 0 hln 100 EDIF=PCB
a 0 u 13 0 0 0 hln 100 SCICARDS=PCB
a 0 u 13 0 0 0 hln 100 POLARIS=PCB
@index
pageloc 1 0 3092 
@status
n 0 97:11:27:22:47:49;946421269 e 
s 0 97:11:27:22:47:58;946421278 e 
*page 1 0 970 720 iA
@ports
@parts
part 1 ABM1 280 250 h
a 0 u 13 0 14 -4 hln 100 EXP1=(V(%IN)/Masse)
a 0 sp 0 0 14 48 hln 100 PART=ABM1
a 0 a 0 0 0 0 hln 100 PKGREF=ABM11
a 0 ap 0 0 -12 -2 hln 100 REFDES=ABM11
part 3 PARAM 80 160 h
a 0 u 13 0 0 20 hln 100 NAME1=Masse
a 0 u 13 0 0 32 hln 100 NAME2=D_ko
a 0 u 13 0 0 46 hln 100 NAME3=F_ko
a 0 u 13 0 92 22 hlb 100 VALUE1=1
a 0 u 13 0 92 34 hlb 100 VALUE2=1.265
a 0 u 13 0 94 46 hlb 100 VALUE3=40
a 0 a 0 0 0 0 hln 100 PKGREF=PM1
a 1 ap 0 0 10 -2 hcn 100 refdes=PM1
part 4 INTEG 390 250 h
a 0 sp 0 0 0 30 hln 100 PART=INTEG
a 0 a 0 0 0 0 hln 100 PKGREF=INTEG1
a 0 ap 0 0 -6 -6 hln 100 REFDES=INTEG1
a 0 u 13 0 20 24 hln 100 GAIN=1.0
part 5 INTEG 470 250 h
a 0 sp 0 0 0 30 hln 100 PART=INTEG
a 0 a 0 0 0 0 hln 100 PKGREF=INTEG2
a 0 ap 0 0 -6 -6 hln 100 REFDES=INTEG2
part 6 DIFF 230 250 h
a 0 sp 0 0 16 32 hln 100 PART=DIFF
a 0 a 0 0 0 0 hln 100 PKGREF=DIFF1
a 1 ap 0 0 10 2 hln 100 REFDES=DIFF1
part 7 DIFF 170 250 h
a 0 sp 0 0 16 32 hln 100 PART=DIFF
a 0 a 0 0 0 0 hln 100 PKGREF=DIFF2
a 1 ap 0 0 10 2 hln 100 REFDES=DIFF2
part 54 ABM1 370 130 H
a 0 u 13 0 14 -6 hln 100 EXP1=(V(%IN) * F_ko)
a 0 sp 0 0 14 48 hln 100 PART=ABM1
a 0 a 0 0 0 0 hln 100 PKGREF=ABM15
a 0 ap 0 0 -12 -2 hln 100 REFDES=ABM15
part 55 ABM1 370 190 H
a 0 u 13 0 14 -6 hln 100 EXP1=(V(%IN) * D_ko)
a 0 sp 0 0 14 48 hln 100 PART=ABM1
a 0 a 0 0 0 0 hln 100 PKGREF=ABM16
a 0 ap 0 0 -12 -2 hln 100 REFDES=ABM16
part 57 ABM 140 250 h
a 0 u 13 0 -4 54 hln 100 EXP1=10*(1-pwr(2.718,-2*Time))
a 0 sp 0 0 14 48 hln 100 PART=ABM
a 0 a 0 0 0 0 hln 100 PKGREF=ABM1
a 0 ap 0 0 -12 -2 hln 100 REFDES=ABM1
part 0 titleblk 970 720 h
a 1 s 13 0 350 10 hcn 100 PageSize=A
a 1 s 13 0 180 60 hcn 100 PageTitle=
a 1 s 13 0 300 95 hrn 100 PageNO=1
a 1 s 13 0 340 95 hrn 100 PageCount=1
part 50 nodeMarker 570 250 h
a 0 a 0 0 4 22 hlb 100 LABEL=1
@conn
w 94
s 140 250 170 250 93
w 96
s 200 250 230 250 95
w 98
s 260 250 280 250 97
w 100
s 370 250 390 250 131
w 102
s 440 250 450 250 101
s 450 250 470 250 111
s 450 250 450 190 109
s 450 190 370 190 112
w 104
s 520 250 570 250 103
s 570 250 570 130 105
s 570 130 370 130 107
w 115
s 280 190 240 190 114
s 240 190 240 240 116
w 119
s 280 130 180 130 118
s 180 130 180 240 120
@junction
j 170 250
+ p 7 IN2
+ w 94
j 140 250
+ p 57 OUT
+ w 94
j 230 250
+ p 6 IN2
+ w 96
j 200 250
+ p 7 OUT
+ w 96
j 280 250
+ p 1 IN
+ w 98
j 260 250
+ p 6 OUT
+ w 98
j 370 250
+ p 1 OUT
+ w 100
j 390 250
+ p 4 in
+ w 100
j 440 250
+ p 4 out
+ w 102
j 470 250
+ p 5 in
+ w 102
j 520 250
+ p 5 out
+ w 104
j 370 130
+ p 54 IN
+ w 104
j 450 250
+ w 102
+ w 102
j 370 190
+ p 55 IN
+ w 102
j 280 190
+ p 55 OUT
+ w 115
j 240 240
+ p 6 IN1
+ w 115
j 280 130
+ p 54 OUT
+ w 119
j 180 240
+ p 7 IN1
+ w 119
j 570 250
+ p 50 pin1
+ w 104
@attributes
a 0 s 0 0 0 0 hln 100 PAGETITLE=
a 0 s 0 0 0 0 hln 100 PAGENO=1
a 0 s 0 0 0 0 hln 100 PAGESIZE=A
a 0 s 0 0 0 0 hln 100 PAGECOUNT=1
@text
s 5 370 280 hln 100 acel
s 5 450 280 hln 100 vel
s 5 540 280 hln 100 disp
s 5 300 250 hln 100 1 / Masse
s 5 300 190 hln 100 Dämpfung
s 5 300 130 hln 100 Feder
s 5 80 260 hln 100 Last
