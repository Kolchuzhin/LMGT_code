*version 6.0 910981200
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
u 876
blk? 27
PM? 5
I? 2
R? 19
@libraries
@analysis
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
pageloc 1 0 7819 
@status
n 0 97:11:27:22:37:26;946420646 e 
s 0 97:11:27:22:37:30;946420650 e 
*page 1 0 970 720 iA
@ports
port 6 AGND 100 450 h
@parts
part 5 Quelle-sta 700 160 u
a 0 a 0 0 0 0 hln 100 PKGREF=I1
a 1 ap 0 0 20 10 hcn 100 REFDES=I1
a 0 u 13 13 -14 20 hlb 100 Kraft={F_z}
part 189 B-Param 110 490 h
a 0 a 0 0 0 0 hln 100 PKGREF=PM3
a 1 ap 0 0 10 -2 hcn 100 refdes=PM3
a 0 u 13 0 58 50 hlb 100 VALUE7=52e9
a 0 ur 13 0 56 30 hlb 100 VALUE1=2329
a 0 u 13 0 52 40 hlb 100 VALUE2=169e9
a 0 u 13 0 48 90 hlb 100 VALUE3=1.07e-9
a 0 u 13 0 44 60 hlb 100 VALUE4=1000e-6
a 0 u 13 0 48 70 hlb 100 VALUE5=8.9e-21
a 0 u 13 0 42 80 hlb 100 VALUE6=3.23e-20
part 533 Param 230 490 h
a 0 u 13 0 0 20 hln 100 NAME1=l_BK
a 0 u 13 0 50 22 hlb 100 VALUE1=333e-6
a 0 a 0 0 0 0 hln 100 PKGREF=PM4
a 1 ap 0 0 10 -2 hcn 100 refdes=PM4
a 0 u 13 0 0 30 hln 100 NAME2=F_z
a 0 u 13 0 48 32 hlb 100 VALUE2=3.75e-5
part 728 2d-beam_BT 130 240 h
a 0 u 13 13 58 52 hlb 100 l={L_BK}
a 0 sp 13 0 22 38 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk19
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk19
part 729 2d-beam_BT 300 290 d
a 0 u 13 13 24 -6 hlb 100 Ex={E_Mod}
a 0 u 13 13 38 -6 hlb 100 Ixx={F_Tr}
a 0 u 13 13 66 -6 hlb 100 l={L_B}
a 0 u 13 13 52 -6 hlb 100 GM={G_Mod}
a 0 u 13 13 80 -6 hlb 100 IT={T_Tr}
a 0 sp 13 0 36 90 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk20
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk20
part 730 2d-beam_BT 350 390 h
a 0 u 13 13 58 52 hlb 100 l={L_BK}
a 0 sp 13 0 22 38 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk21
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk21
part 731 2d-beam_BT 300 120 d
a 0 u 13 13 48 -10 hlb 100 Ex={E_Mod}
a 0 u 13 13 76 -12 hlb 100 Ixx={F_Tr}
a 0 u 13 13 62 -12 hlb 100 l={L_B}
a 0 u 13 13 88 -10 hlb 100 GM={G_Mod}
a 0 sp 13 0 30 98 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk22
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk22
part 732 2d-beam_BT 350 90 h
a 0 u 13 13 58 52 hlb 100 l={L_BK}
a 0 sp 13 0 22 38 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk23
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk23
a 0 u 13 13 58 66 hlb 100 IT={T_Tr}
part 779 2d-beam_BT 520 290 d
a 0 sp 13 0 32 94 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk24
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk24
a 0 u 13 13 62 -6 hlb 100 Ixx={F_Tr}
a 0 u 13 13 48 -6 hlb 100 l={L_B}
a 0 u 13 13 76 -8 hlb 100 GM={G_Mod}
a 0 u 13 13 36 -4 hlb 100 IT={T_Tr}
part 780 2d-beam_BT 520 120 d
a 0 u 13 13 64 -6 hlb 100 Ixx={F_Tr}
a 0 u 13 13 38 -4 hlb 100 l={L_B}
a 0 u 13 13 78 -6 hlb 100 GM={G_Mod}
a 0 u 13 13 52 -4 hlb 100 IT={T_Tr}
a 0 sp 13 0 46 88 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk25
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk25
part 811 2d-beam_BT 580 240 h
a 0 u 13 13 58 52 hlb 100 l={L_BK}
a 0 sp 13 0 22 38 hln 100 PART=2d-beam_BT
a 0 a 0 0 0 0 hln 100 PKGREF=blk26
a 0 ap 0 0 -12 -2 hln 100 REFDES=blk26
part 0 titleblk 970 720 h
a 1 s 13 0 350 10 hcn 100 PageSize=A
a 1 s 13 0 180 60 hcn 100 PageTitle=
a 1 s 13 0 300 95 hrn 100 PageNO=1
a 1 s 13 0 340 95 hrn 100 PageCount=1
part 875 VIEWPOINT 700 240 h
a 0 s 0 0 0 0 hln 100 MARKERTYPE=VIEWPOINT
a 1 s 13 0 14 8 hcn 100 BIASVALUE=2.587E-06
@conn
w 683
a 0 s 0 0 260 238 hln 100 LABEL=n2
s 220 240 300 240 682
a 0 s 3 0 250 236 hln 100 LABEL=n2
s 300 240 300 210 718
s 300 240 300 290 737
w 687
a 0 s 0 0 250 248 hln 100 LABEL=n3
s 220 250 280 250 686
a 0 s 3 0 250 248 hln 100 LABEL=n3
s 280 250 280 210 722
s 280 250 280 290 733
w 691
a 0 s 0 0 255 258 hln 100 LABEL=n4
s 220 260 290 260 690
a 0 s 3 0 249 258 hln 100 LABEL=n4
s 290 260 290 210 720
s 290 260 290 290 735
w 740
a 0 s 0 0 325 388 hln 100 LABEL=n8
s 350 390 300 390 739
a 0 s 3 0 325 388 hln 100 LABEL=n8
s 300 390 300 380 741
w 744
a 0 s 0 0 315 398 hln 100 LABEL=n9
s 350 400 280 400 743
a 0 s 3 0 315 398 hln 100 LABEL=n9
s 280 400 280 380 745
w 748
a 0 s 0 0 320 408 hln 100 LABEL=n10
s 350 410 290 410 747
a 0 s 3 0 320 408 hln 100 LABEL=n10
s 290 410 290 380 749
w 752
a 0 s 0 0 325 88 hln 100 LABEL=n5
s 350 90 300 90 751
a 0 s 3 0 325 88 hln 100 LABEL=n5
s 300 90 300 120 753
w 756
a 0 s 0 0 315 98 hln 100 LABEL=n6
s 350 100 280 100 755
a 0 s 3 0 315 98 hln 100 LABEL=n6
s 280 100 280 120 757
w 760
a 0 s 0 0 320 108 hln 100 LABEL=n7
s 350 110 290 110 759
a 0 s 3 0 320 108 hln 100 LABEL=n7
s 290 110 290 120 761
w 782
a 0 s 0 0 480 388 hln 100 LABEL=n17
s 440 390 520 390 781
a 0 s 3 0 480 388 hln 100 LABEL=n17
s 520 390 520 380 783
w 786
a 0 s 0 0 470 398 hln 100 LABEL=n18
s 440 400 500 400 785
a 0 s 3 0 470 398 hln 100 LABEL=n18
s 500 400 500 380 787
w 790
a 0 s 0 0 475 408 hln 100 LABEL=n19
s 440 410 510 410 789
a 0 s 3 0 475 408 hln 100 LABEL=n19
s 510 410 510 380 791
w 794
a 0 s 0 0 480 88 hln 100 LABEL=n11
s 440 90 520 90 793
a 0 s 3 0 480 88 hln 100 LABEL=n11
s 520 90 520 120 795
w 798
a 0 s 0 0 470 98 hln 100 LABEL=n12
s 440 100 500 100 797
a 0 s 3 0 470 98 hln 100 LABEL=n12
s 500 100 500 120 799
w 802
a 0 s 0 0 475 108 hln 100 LABEL=n13
s 440 110 510 110 801
a 0 s 3 0 475 108 hln 100 LABEL=n13
s 510 110 510 120 803
w 806
a 0 s 0 0 550 238 hln 100 LABEL=n14
s 580 240 520 240 812
a 0 s 3 0 550 238 hln 100 LABEL=n14
s 520 210 520 240 805
s 520 240 520 290 814
w 808
a 0 s 0 0 545 258 hln 100 LABEL=n16
s 580 260 510 260 867
a 0 s 3 0 545 258 hln 100 LABEL=n16
s 510 210 510 260 807
s 510 260 510 290 869
w 810
a 0 s 0 0 540 248 hln 100 LABEL=n15
s 580 250 500 250 864
a 0 s 3 0 540 248 hln 100 LABEL=n15
s 500 210 500 250 809
s 500 250 500 290 866
w 854
a 0 s 3 0 680 98 hln 100 LABEL=n20
s 700 240 670 240 862
a 0 s 3 0 683 198 hln 100 LABEL=n20
s 700 160 700 240 873
w 855
a 0 s 3 0 115 238 hln 100 LABEL=n1
s 100 240 130 240 572
a 0 s 3 0 115 238 hln 100 LABEL=n1
s 700 450 700 260 827
s 670 260 700 260 831
s 100 450 700 450 825
s 100 450 100 260 763
s 100 260 100 250 565
s 130 260 100 260 562
s 100 240 100 250 574
s 130 250 100 250 567
s 100 240 100 50 765
s 700 50 100 50 837
s 670 250 700 250 856
s 700 250 700 260 858
s 700 50 700 120 871
@junction
j 100 250
+ w 855
+ w 855
j 100 240
+ w 855
+ w 855
j 300 240
+ w 683
+ w 683
j 290 260
+ w 691
+ w 691
j 280 250
+ w 687
+ w 687
j 220 240
+ p 728 t_r
+ w 683
j 220 250
+ p 728 r_r
+ w 687
j 220 260
+ p 728 d_r
+ w 691
j 130 250
+ p 728 r_l
+ w 855
j 130 240
+ p 728 t_l
+ w 855
j 130 260
+ p 728 d_l
+ w 855
j 300 210
+ p 731 t_r
+ w 683
j 280 210
+ p 731 d_r
+ w 687
j 290 210
+ p 731 r_r
+ w 691
j 280 290
+ p 729 d_l
+ w 687
j 290 290
+ p 729 r_l
+ w 691
j 300 290
+ p 729 t_l
+ w 683
j 350 390
+ p 730 t_l
+ w 740
j 300 380
+ p 729 t_r
+ w 740
j 350 400
+ p 730 r_l
+ w 744
j 280 380
+ p 729 d_r
+ w 744
j 350 410
+ p 730 d_l
+ w 748
j 290 380
+ p 729 r_r
+ w 748
j 350 90
+ p 732 t_l
+ w 752
j 300 120
+ p 731 t_l
+ w 752
j 350 100
+ p 732 r_l
+ w 756
j 280 120
+ p 731 d_l
+ w 756
j 350 110
+ p 732 d_l
+ w 760
j 290 120
+ p 731 r_l
+ w 760
j 100 450
+ s 6
+ w 855
j 100 260
+ w 855
+ w 855
j 440 390
+ p 730 t_r
+ w 782
j 520 380
+ p 779 t_r
+ w 782
j 440 400
+ p 730 r_r
+ w 786
j 500 380
+ p 779 d_r
+ w 786
j 440 410
+ p 730 d_r
+ w 790
j 510 380
+ p 779 r_r
+ w 790
j 440 90
+ p 732 t_r
+ w 794
j 520 120
+ p 780 t_l
+ w 794
j 440 100
+ p 732 r_r
+ w 798
j 500 120
+ p 780 d_l
+ w 798
j 440 110
+ p 732 d_r
+ w 802
j 510 120
+ p 780 r_l
+ w 802
j 520 290
+ p 779 t_l
+ w 806
j 520 210
+ p 780 t_r
+ w 806
j 510 290
+ p 779 r_l
+ w 808
j 510 210
+ p 780 r_r
+ w 808
j 500 290
+ p 779 d_l
+ w 810
j 500 210
+ p 780 d_r
+ w 810
j 580 240
+ p 811 t_l
+ w 806
j 520 240
+ w 806
+ w 806
j 670 260
+ p 811 d_r
+ w 855
j 670 250
+ p 811 r_r
+ w 855
j 700 260
+ w 855
+ w 855
j 670 240
+ p 811 t_r
+ w 854
j 580 250
+ p 811 r_l
+ w 810
j 500 250
+ w 810
+ w 810
j 580 260
+ p 811 d_l
+ w 808
j 510 260
+ w 808
+ w 808
j 700 120
+ p 5 -
+ w 855
j 700 160
+ p 5 +
+ w 854
j 700 240
+ p 875 +
+ w 854
@attributes
a 0 s 0 0 0 0 hln 100 PAGETITLE=
a 0 s 0 0 0 0 hln 100 PAGENO=1
a 0 s 0 0 0 0 hln 100 PAGESIZE=A
a 0 s 0 0 0 0 hln 100 PAGECOUNT=1
@text
