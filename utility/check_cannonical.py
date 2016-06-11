#! /usr/bin/env python

# splicingDonnorMotif = ["AG", "GTRAGT"]
# splicingAcceptorMotif = ["YYYYNCAG", "G"]

import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        motif = F[19].split(',')
        donnor_acceptor = F[20]
        donnor_acceptor = donnor_acceptor.replace(" disruption", "")
        donnor_acceptor = donnor_acceptor.replace(" creation", "")

        cannonical_start = 0
        cannonical_end = 0

        if motif[1] == '+'

Miya14  3       113321997       113322826       splice-site-slip        SIDT1(NM_017699);SIDT1(NM_001308350)    11;11   e;e     SIDT1(NM_017699);SIDT1(NM_001308350)    12;12   *;*     3       113322805       .       A       G       60      PASS    SOMATIC 3:113322799-113322807,+ splicing acceptor disruption
Miya17  15      63569935        63578757        exon-skip       APH1B(NM_031301);APH1B(NM_001145646)    0;0     e;e     APH1B(NM_031301);APH1B(NM_001145646)    2;2     s;s     15      63571358        .       A       T       60      PASS    SOMATIC 15:63571352-63571360,+  splicing acceptor disruption
Miya2   X       47045029        47045127        splice-site-slip        RBM10(NM_152856);RBM10(NM_001204466);RBM10(NM_001204467);RBM10(NM_001204468);RBM10(NM_005676)   18;18;19;19;19  e;e;e;e;e       RBM10(NM_152856);RBM10(NM_001204466);RBM10(NM_001204467);RBM10(NM_001204468);RBM10(NM_005676)   19;19;20;20;20  *;*;*;*;*       X       47045114        .       G       A       60      PASS    SOMATIC X:47045107-47045115,+   splicing acceptor disruption
Miya23  1       85029101        85031524        exon-skip       CTBS(NM_004388) 5       s       CTBS(NM_004388) 3       e       1       85029416        .       A       G       60      PASS    SOMATIC 1:85029412-85029419,-   splicing donnor disruption
Miya25  1       11150725        11151550        exon-skip       EXOSC10(NM_001001998);EXOSC10(NM_002685)        5;5     s;s     EXOSC10(NM_001001998);EXOSC10(NM_002685)        3;3     e;e     1       11151071        .       C       G       60      PASS    SOMATIC 1:11151065-11151072,-
   splicing donnor disruption
Miya25  3       9516655 9516752 pseudo-exon-inclusion   SETD5(NM_001292043);SETD5(NM_001080517) *;*     *;*     SETD5(NM_001292043);SETD5(NM_001080517) 23;21   s;s     3       9516654 .       G       A       60      PASS    SOMATIC 3:9516654-9516661,+     splicing donnor creation
Miya25  4       5731117 5733221 splice-site-slip        EVC(NM_001306090);EVC(NM_001306092);EVC(NM_153717)      2;2;2   e;e;e   EVC(NM_001306090);EVC(NM_001306092);EVC(NM_153717)      3;3;3   *;*;*   4       5733216 .       A       C       60      PASS    SOMATIC 4:5733213-5733221,+
     splicing acceptor creation


