from scipy import spatial
from scipy.spatial.qhull import QhullError

IMEC2 = False
H32 = True

if IMEC2:
    # site:channel
    s = {0: 0,
         1: 1,
         2: 2,
         3: 3,
         4: 4,
         5: 5,
         6: 6,
         7: 7,
         8: 8,
         9: 9,
         10: 10,
         11: 11,
         12: 12,
         13: 13,
         14: 14,
         15: 15,
         16: 16,
         17: 17,
         18: 18,
         19: 19,
         20: 20,
         21: 21,
         22: 22,
         23: 23,
         24: 24,
         25: 25,
         26: 26,
         27: 27,
         28: 28,
         29: 29,
         }


    channel_groups = {
        # Shank index.
        0: {   
            # List of channels to keep for spike detection.
            'channels': s.values(),

            # 2D positions of the channels
            # channel: (x,y)
            'geometry': {
                    s[0]: (0, 1050),
                    s[1]: (22, 1050),
                    s[2]: (0, 1072),
                    s[3]: (22, 1072),
                    s[4]: (0, 1094),
                    s[5]: (22, 1094),
                    s[6]: (0, 1116),
                    s[7]: (22, 1116),
                    s[8]: (0, 1138),
                    s[9]: (22, 1138),
                    s[10]: (0, 1160),
                    s[11]: (22, 1160),
                    s[12]: (0, 1182),
                    s[13]: (22, 1182),
                    s[14]: (0, 1204),
                    s[15]: (22, 1204),
                    s[16]: (0, 1226),
                    s[17]: (22, 1226),
                    s[18]: (0, 1248),
                    s[19]: (22, 1248),
                    s[20]: (0, 1270),
                    s[21]: (22, 1270),
                    s[22]: (0, 1292),
                    s[23]: (22, 1292),
                    s[24]: (0, 1314),
                    s[25]: (22, 1314),
                    s[26]: (0, 1336),
                    s[27]: (22, 1336),
                    s[28]: (0, 1358),
                    s[29]: (22, 1358),
            }
        }
    }

if H32: 
    # site:channel
    s = {0: 0,
         1: 1,
         2: 2,
         3: 3,
         4: 4,
         5: 5,
         6: 6,
         7: 7,
         }

    channel_groups = {
        # Shank index.
        0: {   
            # List of channels to keep for spike detection.
            'channels': s.values(),

            # 2D positions of the channels
            # channel: (x,y)
            'geometry': {
                    s[0]: (-24, 1100),
                    s[1]: (18, 1130),
                    s[2]: (-18, 1160),
                    s[3]: (12, 1190),
                    s[4]: (-12, 1220),
                    s[5]: (6, 1250),
                    s[6]: (-6, 1280),
                    s[7]: (0, 1310),
            }
        }
    }


print channel_groups


def get_graph_from_geometry(geometry):

    # let's transform the geometry into lists of channel names and coordinates
    chans,coords = zip(*[(ch,xy) for ch,xy in geometry.iteritems()])

    # we'll perform the triangulation and extract the 
    try:
        tri = spatial.Delaunay(coords)
    except QhullError:
        # oh no! we probably have a linear geometry.
        chans,coords = list(chans),list(coords)
        x,y = zip(*coords)
        # let's add a dummy channel and try again
        coords.append((max(x)+1,max(y)+1))
        tri = spatial.Delaunay(coords)

    # then build the list of edges from the triangulation
    indices, indptr = tri.vertex_neighbor_vertices
    edges = []
    for k in range(indices.shape[0]-1):
        for j in indptr[indices[k]:indices[k+1]]:
            try:
                edges.append((chans[k],chans[j]))
            except IndexError:
                # let's ignore anything connected to the dummy channel
                pass
    return edges


def build_geometries(channel_groups):
    for gr, group in channel_groups.iteritems():
        group['graph'] = get_graph_from_geometry(group['geometry'])
    return channel_groups



channel_groups = build_geometries(channel_groups)
print(channel_groups[0]['graph'])
