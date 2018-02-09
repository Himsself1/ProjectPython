from bokeh.io import output_notebook, show, save
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool
import math, decimal
import numpy as np
import re
import pandas as pd
from bokeh.palettes import Spectral11


def flatten( list_of_lists ):
    '''
    concatenates list of lists as a single list
    '''
    return [item for sublist in list_of_lists for item in sublist]
       
def make_pca_plots( eisodos, flags, out ):
    '''
    Takes as input a numpy array and a groupby object and makes a plot
    that is storred as <out> 
    '''
    print( "this is the imported function" )
    data = dict(
        x = eisodos[:,0],
        y = eisodos[:,1],
        size = [12]*eisodos.shape[0],
        colors = flatten([[Spectral11[i]]*j[0] for i, j in enumerate(flags.values)]),
        )
    source = ColumnDataSource( data )
    bliblikia = ['box_zoom','pan','save','reset','tap','wheel_zoom']
    run = figure(
        plot_width = 1000,
        plot_height = 600,
        tools = bliblikia
        )
    run.circle(
        data['x'], data['y'],
        size = data['size'],
        fill_color = data['colors'], 
        fill_alpha = 0.8,
        )
    run.xaxis.axis_label = 'PC1'
    run.yaxis.axis_label = 'PC2'
    save(run, out)
    return( data )


# def read_my_labels( my_vector ):
# 	'''
# 	Takes a pd.Series with label names and returns grouped labels
# 	'''
# 	kefali = pd.DataFrame(my_vector.str.split("_").tolist(), columns=["population", "index"])
# 	kefali["index"] = pd.to_numeric(kefali['index'])
# 	aa = kefali.groupby("population").count()
# 	return aa


    
# data = pd.read_csv( 'pca_file.tsv', delimiter = '\t', names = ['Ind','x','y'], dtype = {'Ind':'object', 'x':np.float64, 'y':np.float64} )


# print(data.shape)
# print(data['x'][0:10])

# print(np.array(data))

# make_pca_plots( np.array(data)[:,1:], read_my_labels(data['Ind']), 'not.plot.html' )
