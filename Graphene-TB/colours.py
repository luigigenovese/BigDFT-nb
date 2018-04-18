#! /usr/bin/env python

import math
import os
import sys

# returns an array of html colour codes for n evenly spaced colours from dark red to violet
def find_colours(n):

    colours={}
    colours['html']=list()
    colours['indices']=list()
    ind3d = [0] * 3
    for i in range(n):
        ind1d=128+i*1408/n
        if (ind1d<256):
          ind3d[0]=ind1d
          ind3d[1]=0
          ind3d[2]=0
        elif (ind1d<512):
          ind3d[0]=255
          ind3d[1]=ind1d-256
          ind3d[2]=0
        elif (ind1d<768):
          ind3d[0]=768-ind1d-1
          ind3d[1]=255
          ind3d[2]=0
        elif (ind1d<1024):
          ind3d[0]=0
          ind3d[1]=255
          ind3d[2]=ind1d-768
        elif (ind1d<1280):
          ind3d[0]=0
          ind3d[1]=1280-ind1d-1
          ind3d[2]=255
        elif (ind1d<=1536):
          ind3d[0]=ind1d-1280
          ind3d[1]=0
          ind3d[2]=255

        colours['indices'].append(ind3d[:])
        #print ind3d

        html=[0] * 6
        htmls=[0] * 6
        html[0]=ind3d[0]/16
        html[1]=(ind3d[0]%16)
        html[2]=ind3d[1]/16
        html[3]=(ind3d[1]%16)
        html[4]=ind3d[2]/16
        html[5]=(ind3d[2]%16)
        for j in range(6):
            if (html[j]==10):
                htmls[j]='A'
            elif (html[j]==11):
                htmls[j]='B'
            elif (html[j]==12):
                htmls[j]='C'
            elif (html[j]==13):
                htmls[j]='D'
            elif (html[j]==14):
                htmls[j]='E'
            elif (html[j]==15):
                 htmls[j]='F'
            else:
                htmls[j]=html[j]
                #write(htmls[j],'(I1)') html[j]
        #print(htmls)
        #print(" ")
        string='#'+str(htmls[0])+str(htmls[1])+str(htmls[2])+str(htmls[3])+str(htmls[4])+str(htmls[5])
        colours['html'].append(string)
        #print(string)
        #print(colours['html'])
        #print(" ")

    #print(colours['indices'])

    return colours




