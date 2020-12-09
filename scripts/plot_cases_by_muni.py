import sys,os
import numpy as np
import pandas as pd
import json
from math import sqrt
from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection

def column(data,col):
    return [row [col] for row in data]

def load_geojson_to_polygons(gjs_file):
    '''Given a geojson, return a dictionary of polygons keyed by country name.

    Adapted from https://github.com/ebov/space-time/blob/master/Scripts/notebooks/auxiliaries/ebov_data.py#L352
    '''
    polygons = {}
    location_points = {}
    locations = []
    handle=pd.read_json(open(gjs_file,'r')) ## load data
    features=handle['features']

    for loc in features: ## iterate through features (locations)
        try:
            poly = np.array(loc['geometry']['coordinates']) ## get coordinates
            location=loc['properties']['name'] ## standardised location name
            locations.append(location)
        except:
            try:
                poly = np.array(loc['geometry']['geometries'])
                location=loc['properties']['name']
            except:
                print(f"bad: {loc['properties']['name']}")


        polygons[location]=[]
        location_points[location]=[]
        if (loc['geometry']['type']=='MultiPolygon') or (loc['geometry']['type']=='GeometryCollection'): ## multiple parts detected
            try:
                for part in poly:
                    for coords in part:
                        xs=column(coords,0)
                        ys=column(coords,1)
                        location_points[location].append(np.vstack(zip(xs,ys)))
            except:
                print(f"bad2: {loc['properties']['name']}")
        if loc['geometry']['type']=='Polygon': ## location is single part
            for coords in poly:
                xs=column(coords,0)
                ys=column(coords,1)
                location_points[location].append(np.vstack(zip(xs,ys)))

        complete_location=[]
        for part in location_points[location]:
            complete_location.append(Polygon(part,True))

        polygons[location]=complete_location

    return polygons, locations

def fix_loc_name(loc):
    if "(" in loc:
        loc=loc.split("(")[0][:-1]
    return loc

def get_centroids_from_polygons(polygons):
    centroids = {}
    for loc in polygons:
        fix_loc_name(loc)
        try:
            coords = polygons[loc][0].get_xy()
            xs = column(coords,0)
            ys = column(coords,1)
            centroids[loc] = ((sum(xs)/len(xs)),(sum(ys)/len(ys)))
        except:
            pass
            # print(f"couldn't handle {loc}: {polygons[loc]}")
    return centroids

def add_background_map(ax, polygons, locations):
    '''
    '''
    for i,loc in enumerate(locations): ## iterate over locations
        # region=country_region[loc] ## get country
        #
        # regionColor=colors[region.lower()] ## get colour map

        regionColor="lightgrey"

        ax.add_collection(PatchCollection(polygons[loc],facecolor=regionColor,edgecolor='black',lw=1,zorder=-1)) ## polygon colour pale

    return ax

def get_case_count(item):
    if "<" in item["CASES"]:
        C = float(item["CASES"][1:])
    else:
        C = float(item["CASES"])
    return C

if __name__ == '__main__':
    gj_fname = "data/belgium-municipalities.geojson"
    case_count_json = "data/COVID19BE_CASES_MUNI_CUM.json"
    # load geojson to define polygons and locations for the map
    polygons, locations = load_geojson_to_polygons(gj_fname)
    centroids = get_centroids_from_polygons(polygons)

    with open(case_count_json, 'r') as f:
        data = json.load(f)

    provinces = set()

    ### Province colors
    cmap = { 'Limburg' : "midnightblue",
             'Liège' : "red",
             'Namur' : "firebrick",
             'VlaamsBrabant' : "dodgerblue",
             'Hainaut' : "tomato",
             'BrabantWallon' : "lightcoral",
             'Luxembourg' : "maroon",
             'Brussels' : "mediumspringgreen",
             'OostVlaanderen' : "royalblue",
             'Antwerpen' : "blue",
             'WestVlaanderen' : "skyblue" }

    fig, ax = plt.subplots(1,1,dpi=500)
    fig.set_figheight(8.5)
    fig.set_figwidth(11)

    ax = add_background_map(ax, polygons, locations)

    point_size_func = lambda x: 2*sqrt(x)
    Alpha = .6
    # print(centroids)
    for item in data:
        try:
            if fix_loc_name(item["TX_DESCR_NL"]) in centroids.keys():
                n = fix_loc_name(item["TX_DESCR_NL"])
                c = centroids[n]
                C = get_case_count(item)
                ax.scatter( c[0],c[1],
                            point_size_func(C),
                            facecolor=cmap[item["PROVINCE"]],
                            alpha=Alpha,
                            zorder=100000
                            )
            elif fix_loc_name(item["TX_DESCR_FR"]) in centroids.keys():
                n = fix_loc_name(item["TX_DESCR_FR"])
                c = centroids[n]
                C = get_case_count(item)
                ax.scatter( c[0],c[1],
                            point_size_func(C),
                            facecolor=cmap[item["PROVINCE"]],
                            alpha=Alpha,
                            zorder=100000
                            )
            else:
                try:
                    n = f"{fix_loc_name(item['TX_DESCR_NL'])}#{fix_loc_name(item['TX_DESCR_FR'])}"
                    c = centroids[n]
                    C = get_case_count(item)
                    ax.scatter( c[0],c[1],
                                point_size_func(C),
                                facecolor=cmap[item["PROVINCE"]],
                                alpha=Alpha,
                                zorder=100000
                                )
                except:
                    print(f"Couldn't find a location to plot {item['TX_DESCR_FR']}")
        except KeyError:
            pass

    ax.set_xlim(2.5,6.5)
    ax.set_ylim(49.45,51.55)

    plt.savefig('test.png')
