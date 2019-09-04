'''
Created on 11 Jan 2019

@author: thomasgumbricht
'''


import os
from sys import exit
import urllib.request
from html.parser import HTMLParser
#import geoimagine.gis.mj_gis_v80 as mj_gis 
from geoimagine.gdalutilities import GDALstuff
from geoimagine.ancillary import ProcessAncillary
from geoimagine.kartturmain import LayerCommon
#from geoimagine.support.karttur_dt import Today
import subprocess

import numpy as np
import gc

GDALpath = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs'

class AncilLayer(LayerCommon):
    '''Class for sentinel tiles'''
    def __init__(self, composition, locusD, datumD, filepath, FN): 
        """The constructor expects an instance of the composition class."""
        LayerCommon.__init__(self)

        self.comp = composition
        
        self.locus = locusD['locus']
        self.locuspath = locusD['path']

        self.path = filepath
        self.FN = FN

        self.datum = lambda: None
        for key, value in datumD.items():
            setattr(self.datum, key, value)
        if self.datum.acqdate:
            self._SetDOY()
            self._SetAcqdateDOY()
        self._SetPath()
        
    def _SetPath(self):
        """Sets the complete path to sentinel tiles"""
        
        self.FP = os.path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.folder, self.locuspath, self.datum.acqdatestr)
        self.FPN = os.path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING FPN contains space %s' %(self.FPN)
            exit(exitstr)

class ProcessDEM:
    '''class for DEM processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        #DEM processing is restricted to tiles only
        print ('        ProcessDEM',self.process.proc.processid) 
        if self.process.proc.processid == 'SRTMsearchDataPool':
            self._DEMSearchDataPool()
        elif self.process.proc.processid == 'SRTMSearchToDB':
            self._SRTMSearchToDB()
        elif self.process.proc.processid == 'downloadSRTMSingleTile':
            self._SRTMSingleTileini()
        elif self.process.proc.processid == 'downloadSRTMRegion':
            self._SRTMRegionIni()
            
        elif self.process.proc.processid[0:3].lower() in ['tpi','tri','rou']:
            self._LoopOneToMany() 
        elif self.process.proc.processid[0:8].lower() in ['landform']:
            self._LoopManyToOne()
        else:
            self._LoopOneToOne()
            
    def _DrillIntoDataPoolSRTM(self): 
        '''
        '''
        self.serverurl = self.process.params.serverurl
        self.version = self.process.params.version
        self.product = self.process.params.product
        if not len(self.version) == 3:
            exit('The version must be 3 digits, e.g. "003" or "006"')
        if not self.version.isdigit():
            exit('The version must be 3 digits, e.g. "003" or "006"')

        #put the remote search path for the requested dates together
        prodPath ='%s.%s' %(self.product,self.version)
        localPath = os.path.join('/volumes',self.process.srcpath.volume,'DataPoolSRTM',prodPath)
        if not os.path.exists(localPath):
            os.makedirs(localPath)
        sensorurl = 'MEASURES'
        return (sensorurl,prodPath,localPath)
            
    def _DEMSearchDataPool(self):
        '''IMPORTANT the user credentials must be in a hidden file in users home directory called ".netrc"
        '''

        sensorurl,prodPath,localPath = self._DrillIntoDataPoolSRTM()
        cmd ='cd %s;' %(localPath)
        os.system(cmd)
        doneFPN = os.path.join(localPath,'done')

        if not os.path.exists(doneFPN):
            os.makedirs(doneFPN)
        #All SRTM data is available under a single date: 2000.02.11
        dateStr = '2000.02.11'

        #search the datapool

        url = os.path.join(self.serverurl,sensorurl,prodPath,dateStr)
        localFPN = os.path.join(localPath,dateStr)
        if self.process.params.product == 'SRTMGL1':
            DEMANDSSUBFOLDERS
        if self.process.params.product == 'SRTMGL1N':
            DEMANDSSUBFOLDERS
         


        if os.path.exists(localFPN) and not self.process.overwrite:
            return
        #Check if the file is in the "done" subfolder
            
        if os.path.exists(os.path.join(doneFPN,dateStr)) and not self.process.overwrite:
            return
        cmd ='cd %s;' %(localPath)
        cmd ='%(cmd)s /usr/local/bin/wget -L --load-cookies --spider --no-parent ~/.datapoolcookies --save-cookies ~/.datapoolcookies %(url)s' %{'cmd':cmd, 'url':url}

        os.system(cmd)
              
    def _SRTMSearchToDB(self):
        '''Load dotapool holdings to local db
            Does not utilize the layer class but take parameters directly from xml
        '''

        prodPath ='%s.%s' %(self.process.params.product, self.process.params.version)

        localPath = os.path.join('/volumes',self.process.srcpath.volume,'DataPoolSRTM',prodPath)
        
        dateStr = '2000.02.11'

        localFPN = os.path.join(localPath,dateStr)
            
        tarFPN = os.path.join(localPath,'done',dateStr)
        if not os.path.exists(os.path.split(tarFPN)[0]):
            os.makedirs(os.path.split(tarFPN)[0])
        if os.path.exists(localFPN):    
            self._ReadSRTMhtml(self.session,localFPN,tarFPN)
        elif os.path.exists(tarFPN):    
            pass                 
        else:
            print ('MODIS bulk file missing', localFPN)
                
    def _ReadSRTMhtml(self,session,FPN,tarFPN):
        tmpFPN,headL = self._ParseSRTMWgetHTML(FPN)
        print (session.name)
        print (tmpFPN)
        session._LoadSRTMBulkTiles(self.process.params,tmpFPN,headL)
        #move the done file to a subdir called done
        move(FPN,tarFPN)
        
    def _ParseSRTMWgetHTML(self, FPN):
        headL = ['tileid','tilefilename','source','product','version','north','south','east','west','lltile','downloaded','organized','exploded']
            
        tmpFP = os.path.split(FPN)[0]
        tmpFP = os.path.split(tmpFP)[0]
        tmpFP = os.path.join(tmpFP,'tmpcsv')
        if not os.path.exists(tmpFP):
            os.makedirs(tmpFP)
        tmpFPN = os.path.join(tmpFP,'tilelist.csv')
        FPN = 'file://%(fpn)s' %{'fpn':FPN}
        req = urllib.request.Request(FPN)
        with urllib.request.urlopen(req) as response:
            html = response.read()
        parser = MjHTMLParser()
        parser.SetLists(headL,self.process.params.product,self.process.params.version)
        parser.feed(str(html)) 
        WriteCSV(parser.hdfL,tmpFPN)
        return tmpFPN, headL
    
    def _SRTMRegionIni(self):
        '''
        '''
        #Always as script as I cannot get curl to run udner python
        sensorurl,prodPath,localPath = self._DrillIntoDataPoolSRTM()
        shFN = '%(prod)s.sh' %{'prod':self.process.params.product}
        shFP = os.path.join(localPath, 'script')
        if not os.path.exists(shFP):
            os.makedirs(shFP)
        shFPN = os.path.join(shFP,shFN)
        self.scriptF = open(shFPN,'w')
        
        
        
        shFN = '%(prod)s.sh' %{'prod':self.process.params.product}
        shFP = os.path.join(localPath, 'script')
        if not os.path.exists(shFP):
            os.makedirs(shFP)
        shFPN = os.path.join(shFP,shFN)
        self.scriptF = open(shFPN,'w')
        #Search all tiles falling insde the defined region
        #Get the extent of the default region
        #print (self.process.proc.userProj.defregion)
        queryD = {'regionid':self.process.proc.userProj.defregion}
        paramL = ['ullon','lrlon','lrlat','ullat']
        
        paramL = ['ullon','lrlon','lrlat','ullat']
        bbox = self.session._SelectDefRegionExtent(queryD, paramL)

        queryD = {'title': '1degsquare'}
        #queryD['lrlon'] = {'op': '>' , 'val': bbox[0]}
        #queryD['ullon'] = {'op': '<' , 'val': bbox[1]}
        #queryD['ullat'] = {'op': '>' , 'val': bbox[2]}
        #queryD['lrlat'] = {'op': '<' , 'val': bbox[3]}
        queryD['ullon'] = bbox[0]
        queryD['lrlon'] = bbox[1]
        queryD['lrlat'] = bbox[2]
        queryD['ullat'] = bbox[3]
        
        #queryD = {'east':self.process.proc.userProj.defregion}
        #queryD = {'htile':h,'vtile':v,'source':source,'product':product,'version':version}
        #queryD['acqdate'] = {'op': '>=', 'val':startdate}
        #SELECT COUNT(*) FROM ancillary.srtmdptiles WHERE south > -25 and north < 25 and west > 0 and east < 10
        paramL = ['regionid']
        tiles = self.session._Select1degSquareTiles(queryD, paramL)
        for tile in tiles:
            self.process.params.lltile = tile[0]
            self._downloadSRTMSingleTile()
        self.scriptF.close()
        print (shFPN)
        FISK
            

    def _SRTMSingleTileini(self):
        #Always as script as I cannot get curl to run udner python
        sensorurl,prodPath,localPath = self._DrillIntoDataPoolSRTM()
        
        shFN = '%(prod)s.sh' %{'prod':self.process.params.product}
        shFP = os.path.join(localPath, 'script')
        if not os.path.exists(shFP):
            os.makedirs(shFP)
        shFPN = os.path.join(shFP,shFN)
        self.scriptF = open(shFPN,'w')
        self._downloadSRTMSingleTile()
        self.scriptF.close()
        print (shFPN)
        FISK
        
    def _downloadSRTMRegion(self):
        SNULLEBAULLE
                
    def _downloadSRTMSingleTile(self):
        '''Download a single tile position
        '''
        print ('    Downloading SRTM tile',self.process.params.lltile)
        
        #Search the data to download
        queryD = {'downloaded': 'N','lltile':self.process.params.lltile,
                  'product':self.process.params.product,
                  'version':self.process.params.version}
        paramL = ['tileid','tilefilename', 'source','product', 'version', 'lltile']
        tile = self.session._SelectSRTMdatapooltilesOntile(queryD, paramL)
        if tile == None:
            SNULLE
            return
        tileD = dict(zip(paramL,tile))
        sensorurl,prodPath,localPath = self._DrillIntoDataPoolSRTM()
        cmd ='cd %s;' %(localPath)
        os.system(cmd)

        #All SRTM data is available under a single date: 2000.02.11
        dateStr = '2000.02.11'

        url = os.path.join(self.serverurl,sensorurl,prodPath,dateStr,tileD['tilefilename'])
        
        localFPN = os.path.join(localPath,tileD['tilefilename'])

        if self.process.params.product == 'SRTMGL1':
            DEMANDSSUBFOLDERS
        if self.process.params.product == 'SRTMGL1N':
            DEMANDSSUBFOLDERS
            
        if os.path.isfile(localFPN):
            #statusD = {'tileid': tile['tileid'],'column':'downloaded', 'status': 'Y'}
            '''
            self.session._UpdateSRTMTileStatus(statusD)
            print (tile['query'])
            self.session._InsertSRTMtile(tile['query'])
            statusD = {'tileid': tile['tileid'],'column':'downloaded', 'status': 'Y'}
            self.session._UpdateSRTMTileStatus(statusD)
            BALLE
            '''
            self._organizeSRTMSingleTile()
            statusD = {'tileid': tile['tileid'],'column':'downloaded', 'status': 'Y'}
            self.session._UpdateSRTMTileStatus(statusD)
            self.session._InsertSRTMtile(tile['query'])

        
        else:
            home = os.path.expanduser("~")
            cookieFPN = os.path.join(home,'.modisdp_cookies')
            cmd = "curl -n -L -c %(c)s -b %(c)s  %(r)s --output %(l)s;" %{'u':self.process.params.remoteuser, 'c':cookieFPN, 'r':url, 'l':localFPN}
            #cmd = "%(cmd)s mv %(output)s %(dstFPN)s;" %{'cmd':cmd,'output':tile['tempFPN'], 'dstFPN':tile['dstFPN']}

            self.scriptF.write(cmd)
            self.scriptF.write('\n')
            '''
            proc=subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
            print ('proc',proc)
            FITTA
            
            commandL = ["curl", "-n", "-L", "-c", cookieFPN, "-b",  cookieFPN, url, "--output",  localFPN]
            subprocess.call(commandL, shell=True)
            subprocess.Popen(commandL, shell=True, stdout=subprocess.PIPE)

            FJULLE
            self._organizeSRTMSingleTile()
            statusD = {'tileid': tile['tileid'],'column':'downloaded', 'status': 'Y'}
            self.session._UpdateSRTMTileStatus(statusD)
            self.session._InsertSRTMtile(tile['query'])
            '''
       
    def _organizeSRTMSingleTile(self):
        #Search the data to organize
        queryD = {'organized': 'N','lltile':self.process.params.lltile,
                  'product':self.process.params.product,
                  'version':self.process.params.version}
        
        paramL = ['tileid','tilefilename', 'source','product', 'version', 'lltile']
        tile = self.session._SelectSRTMdatapooltilesOntile(queryD, paramL)
        if tile == None:
            SNULLE
            return

        tileD = dict(zip(paramL,tile))
        sensorurl,prodPath,localPath = self._DrillIntoDataPoolSRTM()
        '''
        for locus in self.process.srcLayerD:
            print ('srclocus',locus)
            for datum in self.process.srcLayerD[locus]:
                print ('srcdatum',datum)
                for comp in self.process.srcLayerD[locus][datum]:
                    print ('srccomp',comp)
        '''
        
        for locus in self.process.dstLayerD:
            print ('dstlocus',locus)
            if len(self.process.dstLayerD[locus]) == 0:
                exitstr = 'EXITING, no dates defined in Ancillary.OrganizeAncillary'
                exit(exitstr)
            
            for datum in self.process.dstLayerD[locus]:
                print ('dstdatum',datum)
                if len(self.process.dstLayerD[locus][datum]) == 0:                    
                    exitstr = 'EXITING, no compositions defined in Ancillary.OrganizeAncillary'
                    print (exitstr)
                    BALLE
                    exit(exitstr)

                for comp in self.process.dstLayerD[locus][datum]:
    
                    print (self.process.dstLayerD[locus][datum][comp].FPN)
                    print (self.process.dstLayerD[locus][datum][comp].locus.locus)
                    #Change the locus.locus to the tile (instead of global)
                    self.process.dstLayerD[locus][datum][comp].locus.locus = tileD['lltile']
                    self.process.dstLayerD[locus][datum][comp]._SetPath()
                    print (self.process.dstLayerD[locus][datum][comp].FPN)


        self.process.proc.processid = 'organizeancillary'
        print (self.process.proc.paramsD['regionid'])

        
        #self.process.proc.paramsD['regionid'] = self.process.proc.paramsD['lltile']
        ProcessAncillary(self.process,self.session,True)
        '''
        #Construct the ancillary layer
        datadir = self.process.proc.srcraw.paramsD[comp]['datadir']
        self.FP = os.path.join('/Volumes',self.process.srcpath.volume, datadir)
        for FN in os.listdir(self.FP):
            if FN.endswith(self.process.srcpath.hdrfiletype) and solutionStr in FN:
                if FN.endswith('.gz'):
                    zipFPN = os.path.join(self.FP,FN)
                    dstFPN = os.path.splitext(zipFPN)[0]
                    dstFPNbasic,dstext = os.path.splitext(dstFPN)
                    #This is a gunzip file that must be exploded
                    if not os.path.isfile(dstFPN):
                        print ('    unzipping', zipFPN)
                        print ('    to',dstFPN)

                        zipper.ExplodeGunZip(zipFPN)

                    yyyymmdd =  dstFPN.split('.')[2]
                    #Force the date to represent the first day of the month
                    yyyymmdd = '%(yyyymm)s01' %{'yyyymm':yyyymmdd[0:6]}
                    acqdate = mj_dt.yyyymmddDate(yyyymmdd)

                    #Recreate the compositon
                    key = list(self.process.dstCompD.keys())[0]

                    comp = self.process.dstCompD[key]
                    #Here I reset the filetype, this is dangerous 
                    self.process.proc.srcpathD['hdrfiletype'] = dstext
                    
                    #Set the source file (without extension)
                    self.process.proc.srcraw.paramsD[key]['datafile'] = os.path.split(dstFPNbasic)[1]
                      
                    datumD = {'acqdatestr': yyyymmdd, 'acqdate':acqdate}
                
                    filepath = lambda: None
                    filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = 'tif'
                
                    locusD = {'locus':'global', 'path':'global'}
                    #Create a standard reaster layer
                    bandR = RasterLayer(comp, locusD, datumD, filepath)

                    self.process.dstLayerD['global'] = {}
                    self.process.dstLayerD['global'][yyyymmdd] = {key:bandR}
                    
                    srcrawD = self.process.proc.srcraw.paramsD[key]
                    if self.process.proc.srcpathD['hdrfiletype'][0] == '.':
                        ext = self.process.proc.srcpathD['hdrfiletype']
                    else:
                        ext = '.%s' %(self.process.proc.srcpathD['hdrfiletype'])
                    self.srcFN = '%(fn)s%(e)s' %{'fn':srcrawD['datafile'],'e':ext}            
                    self.srcFP = os.path.join('/Volumes',self.process.proc.srcpathD['volume'], srcrawD['datadir'])
                    self.srcFPN = os.path.join(self.srcFP,self.srcFN)
                    ProcessAncillary(self.process,self.session,True)
    '''   
    
            
    def _LoopManyToOne(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                srcCompD = {}
                dstCompL = []
                #Loop the targets for this layer

                for dstComp in self.process.dstLayerD[locus][datum]:

                    if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                        dstCompL.append(dstComp)
                        print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                    elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                        print ('    Done',self.process.dstLayerD[locus][datum][dstComp].FPN)
                        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                if len(dstCompL) == 0:
                    print ('nothing to create')
                    continue

                dstComp = dstCompL[0]

                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if the file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        compid = self.process.srcLayerD[locus][datum][srcComp].comp.id
                        srcCompD[compid] = srcComp
                        srcCompD[compid] = self.process.srcLayerD[locus][datum][srcComp]
                        srcCompD[compid].ReadRasterLayer()
                    else:
                        print ('missing',self.process.srcLayerD[locus][datum][srcComp].FPN)
                        SNULLEBULLE
                        break
                        
                if self.process.proc.processid[0:8].lower() == 'landform':
                    self._LandformTPIini(locus,datum,srcCompD,dstComp)
                
                else:
                    print (self.process.proc.processid)
                    SNULLE

                #Register the layer
                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)
                self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = None
                self.process.dstLayerD[locus][datum][dstComp] = None
                for comp in srcCompD:
                    srcCompD[comp].layer.NPBAND = None
                    srcCompD[comp] = None
        
                #grabage collect
                gc.collect()
            
    def _LoopOneToMany(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstCompL = []
                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if the file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        #Loop the targets for this layer
                        for dstComp in self.process.dstLayerD[locus][datum]:

                            if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                                dstCompL.append(dstComp)
                                print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                            elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                        if len(dstCompL) == 0:
                            continue
                        #Get the raster metadata for the source layer
                        self.process.srcLayerD[locus][datum][srcComp].GetRastermetadata()
                        if self.process.proc.processid[0:3].lower() == 'tpi':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'tpi')
                        elif self.process.proc.processid[0:3].lower() == 'tri':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'tri')
                        elif self.process.proc.processid[0:9].lower() == 'roughness':
                            self._TPITRI(locus,datum,srcComp,dstCompL,'rn')
                        elif self.process.proc.processid[0:9].lower() == 'hillshade':
                            self._HillShade(locus,datum,srcComp,dstCompL)
                        else:
                            print (self.process.proc.processid)
                            SNULLE
                        for dstcomp in dstCompL:
                            self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)
                      
    def _LoopOneToOne(self):
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstCompL = []
                for srcComp in self.process.srcLayerD[locus][datum]:
                    #Check if src file is there
                    if os.path.exists(self.process.srcLayerD[locus][datum][srcComp].FPN):
                        #Loop the targets for this layer
                        for dstComp in self.process.dstLayerD[locus][datum]:
                            if not self.process.dstLayerD[locus][datum][dstComp]._Exists() or self.process.overwrite:
                                dstCompL.append(dstComp)
                                print ('    Processing',self.process.dstLayerD[locus][datum][dstComp].FPN)
                            elif self.process.dstLayerD[locus][datum][dstComp]._Exists():
                                self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)

                        if len(dstCompL) == 0:
                            continue
                        #Get the raster metadata for the source layer
                        self.process.srcLayerD[locus][datum][srcComp].GetRastermetadata()
                        if self.process.proc.processid[0:9].lower() == 'hillshade':
                            self._HillShade(locus,datum,srcComp,dstCompL)
                        elif self.process.proc.processid[0:5].lower() == 'slope':
                            self._Slope(locus,datum,srcComp,dstCompL)
                        elif self.process.proc.processid[0:6].lower() == 'aspect':
                            self._Aspect(locus,datum,srcComp,dstCompL)
                        else:
                            print (self.process.proc.processid)
                            SNULLE
                        for dstcomp in dstCompL:
                            self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)
                     
    def _TPITRI(self,locus,datum,srcComp,dstCompL,txi):
        lins = self.process.srcLayerD[locus][datum][srcComp].comp.metadata.lins
        cols = self.process.srcLayerD[locus][datum][srcComp].comp.metadata.cols
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility

        for dstComp in dstCompL:
            ot =  self.process.dstLayerD[locus][datum][dstComp].comp.celltype
            #The principal steps include 1. changing resolution, 2. Running analysis, 
            #3 resampling back to original resolutuon, and 4 cut out tile
            if self.process.proc.tpiD[dstComp]['resolfac'] == 1:
                scaledSrcFPN = srcFPN

                dsttempFP = os.path.split(self.process.dstLayerD[locus][datum][dstComp].FPN)[0]
                scaledDstFPN = os.path.join(dsttempFP,'scaledTPI.tif')
            else:
                xsize = int(cols/self.process.proc.tpiD[dstComp]['resolfac'])
                ysize = int(lins/self.process.proc.tpiD[dstComp]['resolfac'])
                srctempFP = os.path.split(srcFPN)[0]
                scaledSrcFPN = os.path.join(srctempFP,'scaledDEM.tif')
                dsttempFP = os.path.split(self.process.dstLayerD[locus][datum][dstComp].FPN)[0]
                scaledDstFPN = os.path.join(dsttempFP,'scaledTPI.tif')
                gdaltranslate = GDALstuff(srcFPN, scaledSrcFPN, self.process.params)
                gdaltranslate.TransformOutSize(xsize,ysize,'average',ot)
                
            gdaltxi = GDALstuff(scaledSrcFPN, scaledDstFPN, self.process.params)
            if txi == 'tpi':
                #Run the GDAL TPI utility     
                gdaltxi.TPI()
            elif txi == 'tri':
                #Run the GDAL TPI utility   
                gdaltxi.TRI()
            elif txi == 'rn':
                #Run the GDAL TPI utility   
                gdaltxi.Roughness()
            else:
                exit('unknown terrain index')
            #Copy dst layer and clean up
            if self.process.proc.tpiD[dstComp]['resolfac'] != 1:
                #Resample back to the original size and delete intermediate layers
                xsize = cols
                ysize = lins
                gdaltranslate = GDALstuff(scaledDstFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
                gdaltranslate.TransformOutSize(xsize,ysize,'nearest',ot)
                os.remove(scaledSrcFPN)
            else:
                gdaltranslate = GDALstuff(scaledDstFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
                gdaltranslate.TransformOT(ot)
            os.remove(scaledDstFPN)
                
    def _HillShade(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalHillShade = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalHillShade.HillShade()
      
    def _Slope(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalSlope = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalSlope.Slope()
            
    def _Aspect(self,locus,datum,srcComp,dstCompL):
        if self.process.params.mosaic:
            srcFPN = self._MosaicAdjacentTiles()
        else:
            srcFPN = self.process.srcLayerD[locus][datum][srcComp].FPN
        #The analysis as such is run as GDAL utility
        for dstComp in dstCompL: 
            gdalAspect = GDALstuff(srcFPN, self.process.dstLayerD[locus][datum][dstComp].FPN, self.process.params)
            gdalAspect.Aspect()
                   
    def _LandformTPIini(self,locus,datum,srcCompD,dstComp):
        '''
        '''
        standardize = self.process.params.standardize
        expanded = self.process.params.expanded
        tpiTH = self.process.params.tpithreshold
        slopeTH = self.process.params.slopethreshold
        tpisstd = self.process.params.tpisstd
        tpilstd = self.process.params.tpilstd
        for comp in srcCompD:
            if comp == 'tpis':
                TPIS = srcCompD[comp].layer.NPBAND
                tpisnull = srcCompD[comp].layer.cellnull
            elif comp =='tpim':
                TPIL = srcCompD[comp].layer.NPBAND
                tpilnull = srcCompD[comp].layer.cellnull
            elif comp =='slope':
                SLOPE = srcCompD[comp].layer.NPBAND
                slopenull = srcCompD[comp].layer.cellnull
            else:
                exit('unknown input band in landform')
        #Create the dst array
        dstArr = np.empty_like(TPIS, dtype=np.uint8)
        
        if standardize:
            #Standardize the TPI prior to Landform classification
            if tpisstd and tpilstd:
                print ('    Global standardized')

                self._LandformGlobalStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
            elif tpisstd and not tpilstd:
                print ('    Relative standardized')

                self._LandformRelativeStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
            else:
                print ('    Local standardized')

                #Calcualte local std for TPI (large and small)
                self._LandformStandardized(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded)
      
        else:
            print ('    Raw TPI classification')

            #Landform classification against original TPI data
            if expanded:
                self._LandformTPI12(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH)
            else:
                #Landform classification against original TPI data
                self._LandformTPI10(TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH)

        dstArr[np.where( TPIS == tpisnull )] = 255
        dstArr[np.where( TPIL == tpilnull )] = 255
        dstArr[np.where( SLOPE == slopenull )] = 255

        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstArr
        
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(srcCompD[comp].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()
        
    def _LandformStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        TPISf = TPIS.astype(np.float32)
        TPISf[TPISf==tpisnull] = np.nan
        
        tpiSmean = np.nanmean(TPISf)
        tpiSstd = np.nanstd(TPISf)
        
        print ('    small TPI mean std',tpiSmean, tpiSstd)
        TPISstandard = ((TPISf-tpiSmean)/tpiSstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)
        TPILf[TPILf==tpilnull] = np.nan
        
        tpiLmean = np.nanmean(TPILf)
        tpiLstd = np.nanstd(TPILf)
        print ('    large TPI mean std',tpiLmean, tpiLstd)

        TPILstandard = ((TPILf-tpiLmean)/tpiLstd)*50+0.5
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
                   
    def _LandformGlobalStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        tpisstd = self.process.params.tpisstd
        tpilstd = self.process.params.tpilstd 
        TPISf = TPIS.astype(np.float32)

        TPISstandard = (TPISf/tpisstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)

        TPILstandard = (TPILf/tpilstd)*50+0.5
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
            
    def _LandformRelativeStandardized(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH,tpisnull,tpilnull,expanded):
        '''
        '''
        tpisstd = self.process.params.tpisstd
        
        TPISf = TPIS.astype(np.float32)
        TPISf[TPISf==tpisnull] = np.nan

        TPISstandard = (TPISf/tpisstd)*50+0.5
        
        TPILf = TPIL.astype(np.float32)
        TPILf[TPILf==tpilnull] = np.nan
        
        
        Sstd = np.nanstd(TPISf)
        
        Lstd = np.nanstd(TPILf)
        relativeSD =  Lstd*tpisstd/Sstd
        
        TPILstandard = (TPILf/relativeSD)*50+0.5
                
        if expanded:
            self._LandformTPI12(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
        else:
            #Landform classification against original TPI data
            self._LandformTPI10(TPISstandard,TPILstandard,SLOPE,dstArr,tpiTH,slopeTH)
            
    def _LandformTPIWeissOriginal(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH):
                
        #Plain = 5
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 50

        #Open slope = 6
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 60
        
        # Mesa or flat ridge = 7
        #upper slope edge = 71
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 71
        #mesa or flat ridge top = 72
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 72
        
        #U-shaped valley = 4        
        #valley floor edge = 41
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) & (SLOPE > slopeTH))] = 41
        #central valley floor  = 42
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) & (SLOPE <= slopeTH))] = 42
        
        # shallow valley / midslope drainage = 2
        #shallow valley edge = 21
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 21
        #central shallow valley = 21
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 22
        
        #midslope ridge / hills in valleys = 9
        #midslope ridge / hills in valleys = 91
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 91
        #midslope ridge / hills in valleys = 92
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 92
        
        #upland drainage = 3
        #upland drainage slope = 31
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 31
        #upland drainage flat = 32
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 32
        
        #canyon / incised stream = 1
        #incised canyon - rapid = 11
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) & (SLOPE > slopeTH))] = 11
        #incised canyon - calm = 12
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) & (SLOPE <= slopeTH))] = 12
        
        #peak / mt top = 10
        #peak / mt top = 101
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) & (SLOPE > slopeTH))] = 101
        #peak / mt top = 102
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) & (SLOPE <= slopeTH))] = 102
        
        #local ridge = 8
        #local ridge = 81
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL < -tpiTH) & (SLOPE > slopeTH))] = 81
        #local ridge = 82
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL < -tpiTH) & (SLOPE <= slopeTH))] = 82  
        
    def _LandformTPI10(self,TPIS,TPIL,SLOPE,dstArr,tpiTH,slopeTH):
        '''
        '''
        #Plain 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH) )] = 10
        #Open slope
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH) )] = 11
        # Mesa or flat ridge
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH) )] = 12
        #U-shaped valley       
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) )] = 9
        # shallow valley / midslope drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) )] = 6
        #midslope ridge / hills in valleys
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) )] = 14
        #upland drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH) )] = 8
        #canyon / incised stream
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) )] = 5
        #peak / mt top
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) )] = 16
        #local ridge
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL <= -tpiTH) )] = 13
        
    def _LandformTPI12(self, TPIS, TPIL, SLOPE, dstArr,tpiTH,slopeTH):
        '''
        '''
        #Incised stream /canyon
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL <= -tpiTH) )] = 5
        #channel / shallow valley
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE < slopeTH))] = 6
        #Midslope drainage
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE >= slopeTH))] = 7
        #Upland drainage, Headwater 
        dstArr[np.where( (TPIS <= -tpiTH) & (TPIL >= tpiTH))] = 8
        #U-shaped valley 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL <= -tpiTH) )] = 9
        #Plain 
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 10
        #Open slope = 6, changed to 7
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE > slopeTH))] = 11
        #mesa or flat ridge top = 8 (should include peat dome edges)
        dstArr[np.where( (TPIS > -tpiTH) & (TPIS < tpiTH) & (TPIL >= tpiTH))] = 12
        #local ridge
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL <= -tpiTH) )] = 13
        #hills in valleys 
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE < slopeTH))] = 14
        #midslope ridge = 92, new = 11
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL > -tpiTH) & (TPIL < tpiTH) & (SLOPE <= slopeTH))] = 15
        #Mountain top
        dstArr[np.where( (TPIS >= tpiTH) & (TPIL >= tpiTH) )] = 16
        
class MjHTMLParser(HTMLParser):
    
    def SetLists(self,headL,product,version):
        self.hdfL = []
        self.xmlL = []   
        self.hdfL.append(headL)
        self.xmlL.append(headL)
        self.product = product
        self.version = version
        
    
    def SplitSRTMFileName(self,value):

        lltile = value.split('.')[0]
        if 'E' in lltile:
            n,e = lltile.split('E')
            e = int(e)
            n = int(n[1:len(n)])
            if 'S' in lltile: 
                n *= -1
        else:
            n,e = lltile.split('W')
            e = int(e)
            n = int(n[1:len(n)])
            e *= -1
            if 'S' in lltile: 
                n *= -1

        east = e+1
        west = e
        north = n+1
        south = n      
        #if 'S' in lltile: 
        #    print (lltile,n,e)
        '''
            tileid character(26),
            tilefilename character varying(128),
            source character varying(32) DEFAULT 'SRTM', 
            product varchar(24), 
            version character(3),
            north smallint,
            south smallint,
            east smallint,
            west smallint,
            lltile char(7),
        '''
        #source = 'SRTM'

        tileid = '%(prod)s-%(v)s-%(hv)s' %{'prod':self.product,'v':self.version, 'hv':lltile }
        hdfL = [tileid, value,'SRTM', self.product, self.version, north, south, east, west, lltile,'N','N','N']
        #D = {'tileid':tileid,'version':version,'tilefilename':value,'source':'MODIS','product':product,'acqdate':acqdate,'doy':doy,'folder':'orignal','htile':htile,'vtile':vtile}
        self.hdfL.append(hdfL)
        
    def handle_starttag(self, tag, attrs):
        # Only parse the 'anchor' tag.
        if tag == "a":
            # Check the list of defined attributes.
            for name, value in attrs:
                # If href is defined, print it.
                if name == "href":
                    ext = os.path.splitext(value)[1]
                    if ext.lower() == '.nc':
                        self.SplitSRTMFileName(value)
                        
def WriteCSV(csvL,tmpFPN):
    import csv
    with open(tmpFPN, 'w') as csvfile:
        wr = csv.writer(csvfile)
        wr.writerows(csvL)
        '''
        for row in enumerate(csvL):
            print (row)

            if x > 0:
                hvtile = row[0].split('-')[3]
                h = int(hvtile[1:3])
                v = int(hvtile[4:6])

            wr.writerow(row)
        '''