'''
Created on 11 Jan 2019

@author: thomasgumbricht
'''
__version__ = '0.3.1'
VERSION = tuple( int(x) for x in __version__.split('.') )
metadataD = { 'name':'dem', 
             'author':'Thomas Gumbricht', 
             'author_email':'thomas.gumbricht@gmail.com',
             'title':'Digital Elevation Model processing.', 
             'label':'Digital Elevation Model processing.',
             'prerequisites':'GDAL must be installed',
             'image':'avg-trmm-3b43v7-precip_3B43_trmm_2001-2016_A'}