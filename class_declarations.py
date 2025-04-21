from warnings import filterwarnings
filterwarnings('ignore')
import osmium.geom
import osmium
from osgeo import ogr
from osgeo import osr
import sys
import geopandas as gpd
import math
import os
import shapely
import networkx as nx
from copy import deepcopy


class OsmiumGeom:

    def __init__(self, path_to_pbf: str, polygons_id: list, buffer=True) -> None:
        self.__pbf_path = path_to_pbf
        distr_geom = []
        self.__wkt_fact = osmium.geom.WKTFactory()
        self.__restrictions = {}

        for rel_proc in osmium.FileProcessor(self.__pbf_path).with_filter(osmium.filter.EntityFilter(osmium.osm.RELATION)).with_locations():
            if rel_proc.id in polygons_id:
                local_border = []
                for mbs in rel_proc.members:
                    if mbs.role == "outer":
                        if mbs.type == "w":
                            outer_linestring = []
                            inner_way = osmium.FileProcessor(self.__pbf_path)\
                                .with_filter(osmium.filter.IdFilter([mbs.ref])).with_locations()
                            for iwy in inner_way:
                                try:
                                    outer_linestring.append(shapely.from_wkt(self.__wkt_fact.create_linestring(iwy)))
                                except RuntimeError:
                                    continue

                            unary_line = shapely.unary_union(outer_linestring)
                            if unary_line.geom_type == "LineString":
                                geometry_line = list(unary_line.coords)
                                for local_facet in geometry_line:
                                    if local_facet not in local_border:
                                        local_border.append(local_facet)
                            elif unary_line.geom_type == "MultiLineString":
                                for child_type in unary_line:
                                    geometry_line = list(child_type.coords)
                                    for local_facet in geometry_line:
                                        if local_facet not in local_border:
                                            local_border.append(local_facet)

                if local_border[0] != local_border[-1]:
                    local_border.append(local_border[0])
                local_boundary = shapely.Polygon(local_border)
                local_boundary = shapely.make_valid(local_boundary)
                distr_geom.append(local_boundary)
            else:
                road_restriction = rel_proc.tags.get("restriction")
                if road_restriction is not None:
                    local_restriction = {}
                    for memb in rel_proc.members:
                        if memb.type == "w":
                            local_restriction[memb.role] = memb.ref
                            local_restriction["restriction"] = road_restriction
                    if len(local_restriction) > 0:
                        if 'from' in local_restriction and 'to' in local_restriction:
                            restriction_nodes = (int(local_restriction['from']), int(local_restriction["to"]))
                            self.__restrictions[restriction_nodes] = road_restriction
                
        if buffer:
            shapely_concave = shapely.to_wkt(shapely.concave_hull(shapely.make_valid(shapely.unary_union(distr_geom)), 0.5))
            srs_wgs = osr.SpatialReference()
            srs_wgs.ImportFromEPSG(4326)
            srs_proj = osr.SpatialReference()
            srs_proj.ImportFromEPSG(32637)
            coord_transform = osr.CoordinateTransformation(srs_wgs, srs_proj)

            ogr_geom : ogr.Geometry = ogr.CreateGeometryFromWkt(shapely_concave)
            ogr_geom.AssignSpatialReference(srs_wgs)
            ogr_geom.SwapXY()
            ogr_geom.TransformTo(srs_proj)
            ogr_geom.Buffer(1000)
            ogr_geom.TransformTo(srs_wgs)
            ogr_geom.SwapXY()
            self.__district_geometry = shapely.from_wkt(ogr_geom.ExportToWkt())
        else:
            self.__district_geometry = shapely.concave_hull(shapely.make_valid(shapely.unary_union(distr_geom)), 0.5)


    def __load_roads(self) -> None:
        counter1 = 0
        counter2 = 0
        self.__geometry_list = []
        bad_highway_tags = ["pedestrian", "footway", "bridleway", "steps", "corridor", "path", "sidewalk", "crossing", "traffic_island", "cycleway"]
        for osm_way in osmium.FileProcessor(self.__pbf_path).with_filter(osmium.filter.EntityFilter(osmium.osm.WAY)).with_locations():
            if osm_way.tags.get("highway"):
                if osm_way.tags.get("highway") not in bad_highway_tags:
                    counter1 += 1
                    if osm_way is not None:
                        try:
                            linestr_shapely: shapely.LineString = shapely.from_wkt(self.__wkt_fact.create_linestring(osm_way))
                        except RuntimeError:
                            pass
                        if (linestr_shapely.is_valid):
                            if (linestr_shapely.intersects(self.__district_geometry)):
                                single_road = {}
                                single_road["osm_id"] = osm_way.id
                                single_road["highway"] = osm_way.tags.get("highway")
                                single_road["oneway"] = osm_way.tags.get("oneway")
                                single_road["geometry"] = linestr_shapely
                                self.__geometry_list.append(single_road)
                                counter2 += 1


    
    def send_restrictions(self) -> dict:
        return self.__restrictions
                

    @staticmethod
    def haversine_distance(cds_1: list, cds_2: list) -> float:

        lat1 = cds_1[0]
        lat2 = cds_2[0]
        lon1 = cds_1[1]
        lon2 = cds_2[1]

        lat_diff = (lat2 - lat1) * math.pi / 180
        lat_semi = ((lat1 + lat2) / 2) * math.pi / 180
        lon_diff = (lon2 - lon1) * math.pi / 180
        earth_radious = 6371.0

        first_sin = pow(math.sin(lat_diff / 2), 2)
        second_sin = 1 - pow(math.sin(lat_diff / 2), 2) - pow(math.sin(lat_semi), 2)
        third_sin = pow(math.sin(lon_diff / 2), 2)
        return 1000 * 2 * earth_radious * math.asin(math.sqrt(first_sin + (second_sin * third_sin)))


    def __split_osm(self) -> None:

        self.__splitted_list_view = []

        for val in self.__geometry_list:
            linear_geom = list(val["geometry"].coords)
            for i in range(len(linear_geom) - 1):
                inner_dict = {}

                first_wgs = [round(linear_geom[i][0], 5), round(linear_geom[i][1], 5)]
                second_wgs = [round(linear_geom[i+1][0], 5), round(linear_geom[i+1][1], 5)]
                coord_seq_wgs = shapely.LineString([first_wgs, second_wgs])

                inner_dict['osm_id'] = val['osm_id']
                inner_dict['highway'] = val['highway']
                inner_dict['oneway'] = val['oneway']
                inner_dict['meter_dist'] = round(OsmiumGeom.haversine_distance(first_wgs, second_wgs), 2)
                inner_dict['geometry'] = coord_seq_wgs
                self.__splitted_list_view.append(inner_dict)
                


    def __index_calc(self) -> None:
        self.__node_dict = {}
        node_index = 0
        processed_nodes = set()
        for val in self.__splitted_list_view:
            for tupvalue in list(val["geometry"].coords):
                if (tupvalue not in processed_nodes):
                    self.__node_dict[tupvalue] = node_index
                    processed_nodes.add(tupvalue)
                    node_index += 1

    def __prepare_graph(self) -> list:
        to_nx_list = []

        for val in self.__splitted_list_view:
            u = self.__node_dict[list(val["geometry"].coords)[0]]
            v = self.__node_dict[list(val["geometry"].coords)[1]]
            if (val["oneway"] == "no") or (val["oneway"] == "None") or (val["oneway"] is None):

                nested_dict1 = {}
                nested_dict1["u"] = u
                nested_dict1["v"] = v
                nested_dict1["osm_id"] = val["osm_id"]
                nested_dict1["highway"] = val["highway"]
                nested_dict1["oneway"] = val["oneway"]
                nested_dict1["meter_dist"] = val["meter_dist"]
                nested_dict1["geometry"] = val["geometry"]
                to_nx_list.append(nested_dict1)

                nested_dict2 = {}
                nested_dict2["u"] = v
                nested_dict2["v"] = u
                nested_dict2["osm_id"] = val["osm_id"]
                nested_dict2["highway"] = val["highway"]
                nested_dict2["oneway"] = val["oneway"]
                nested_dict2["meter_dist"] = val["meter_dist"]
                nested_dict2["geometry"] = val["geometry"]
                to_nx_list.append(nested_dict2)

            elif (val["oneway"] == "yes"):

                nested_dict1 = {}
                nested_dict1["u"] = u
                nested_dict1["v"] = v
                nested_dict1["osm_id"] = val["osm_id"]
                nested_dict1["highway"] = val["highway"]
                nested_dict1["oneway"] = val["oneway"]
                nested_dict1["meter_dist"] = val["meter_dist"]
                nested_dict1["geometry"] = val["geometry"]
                to_nx_list.append(nested_dict1)

            elif ((val["oneway"] == "-1") or (val["oneway"] == -1)):
                nested_dict2 = {}
                nested_dict2["u"] = v
                nested_dict2["v"] = u
                nested_dict2["osm_id"] = val["osm_id"]
                nested_dict2["highway"] = val["highway"]
                nested_dict2["oneway"] = val["oneway"]
                nested_dict2["meter_dist"] = val["meter_dist"]
                nested_dict2["geometry"] = val["geometry"]
                to_nx_list.append(nested_dict2)
        
        return to_nx_list
    
    def prepare_osm(self) -> list:
        self.__load_roads()
        self.__split_osm()
        self.__index_calc()
        return self.__prepare_graph()


class nxGraph:
    def __init__(self, prepared_data: list) -> None: # Должен быть в виде geojson
        self.__geodata = prepared_data
    
    def get_restrictions(self, restrct: dict) -> None:
        self.__restrictions = restrct
    
    def create_graph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        for elem in self.__geodata:
            g.add_edge(elem["u"], elem["v"])
        return g

    def __make_dual_noturns(self) -> None:

        def append_to_dict(first_key, first_value, second_key, second_value) -> dict:
            nxsingle_node = {}
            nxsingle_node["u"] = first_key
            nxsingle_node["v"] = second_key
            nxsingle_node["meters"] = first_value["meter_dist"] / 2 + second_value["meter_dist"] / 2
            return nxsingle_node

        edge_id = 0
        self.__reverse_data_noturns = {}
        self.__nxdual_noturns = []
        processed_data = set()
        for ge in self.__geodata:
            proc_nodes = (ge["u"], ge["v"])
            if proc_nodes not in processed_data:
                if (ge["oneway"] == "no") or (ge["oneway"] == "None") or (ge["oneway"] is None):
                    proc_nodes1 = (ge["u"], ge["v"])
                    processed_data.add(proc_nodes1)
                    proc_nodes2 = (ge["v"], ge["u"])
                    processed_data.add(proc_nodes2)
                    self.__reverse_data_noturns[edge_id] = ge
                    edge_id += 1
                else:
                    proc_nodes1 = (ge["u"], ge["v"])
                    processed_data.add(proc_nodes1)
                    self.__reverse_data_noturns[edge_id] = ge
                    edge_id += 1
        
        for gkey, gval in self.__reverse_data_noturns.items():
            for nkey, nval in self.__reverse_data_noturns.items():
                if (gval["oneway"] == "no") or (gval["oneway"] == "None") or (gval["oneway"] is None):
                    if gval["v"] == nval["u"]:
                        self.__nxdual_noturns.append(append_to_dict(gkey, gval, nkey, nval))
                    if gval["u"] == nval["v"]:
                        self.__nxdual_noturns.append(append_to_dict(gkey, gval, nkey, nval))
                        
                    
                else:
                    if gval["v"] == nval["u"]:
                        self.__nxdual_noturns.append(append_to_dict(gkey, gval, nkey, nval))
                        

    def __make_dual_withrturns(self) -> None:

        def append_to_dict(first_key, first_value, second_key, second_value) -> dict:
            nxsingle_node = {}
            nxsingle_node["u"] = first_key
            nxsingle_node["v"] = second_key
            nxsingle_node["meters"] = first_value["meter_dist"] / 2 + second_value["meter_dist"] / 2
            return nxsingle_node

        include_num = 0
        edge_id = 0
        self.__reverse_data_rturns = {}
        self.__nxdual_rturns = []
        processed_data = set()
        for ge in self.__geodata:
            proc_nodes = (ge["u"], ge["v"])
            if proc_nodes not in processed_data:
                if (ge["oneway"] == "no") or (ge["oneway"] == "None") or (ge["oneway"] is None):
                    proc_nodes1 = (ge["u"], ge["v"])
                    processed_data.add(proc_nodes1)
                    proc_nodes2 = (ge["v"], ge["u"])
                    processed_data.add(proc_nodes2)
                    self.__reverse_data_rturns[edge_id] = ge
                    edge_id += 1
                else:
                    proc_nodes1 = (ge["u"], ge["v"])
                    processed_data.add(proc_nodes1)
                    self.__reverse_data_rturns[edge_id] = ge
                    edge_id += 1
        
        for gkey, gval in self.__reverse_data_rturns.items():
            for nkey, nval in self.__reverse_data_rturns.items():
                if (gval["oneway"] == "no") or (gval["oneway"] == "None") or (gval["oneway"] is None):
                    graph_tuple1 = (int(gval["osm_id"]), int(nval["osm_id"]))
                    graph_tuple2 = (int(gval["osm_id"]), int(nval["osm_id"]))
                    if graph_tuple1 in self.__restrictions or graph_tuple2 in self.__restrictions:
                        try:
                            if (self.__restrictions[graph_tuple1].find("no_") != -1) and (self.__restrictions[graph_tuple1].find("u_turn") == -1):
                                include_num += 1
                            elif (self.__restrictions[graph_tuple1].find("no_") != -1) and (self.__restrictions[graph_tuple1].find("u_turn") != -1):
                                include_num += 1
                                if gval["v"] == nval["u"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                            elif self.__restrictions[graph_tuple1].find("only_") != -1:
                                include_num += 1
                                if gval["v"] == nval["u"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                                if gval["u"] == nval["v"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                                else:
                                    pass
                        except KeyError:
                            if (self.__restrictions[graph_tuple2].find("no_") != -1) and (self.__restrictions[graph_tuple2].find("u_turn") == -1):
                                include_num += 1
                            elif (self.__restrictions[graph_tuple2].find("no_") != -1) and (self.__restrictions[graph_tuple2].find("u_turn") != -1):
                                include_num += 1
                                if gval["v"] == nval["u"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                            elif self.__restrictions[graph_tuple2].find("only_") != -1:
                                include_num += 1
                                if gval["v"] == nval["u"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                                if gval["u"] == nval["v"]:
                                    self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                                else:
                                    pass
                    else:
                        if gval["v"] == nval["u"]:
                            self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                        if gval["u"] == nval["v"]:
                            self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                       
                else:
                    graph_tuple = (int(gval["osm_id"]), int(nval["osm_id"]))
                    if graph_tuple in self.__restrictions:
                        if (self.__restrictions[graph_tuple].find("no_") != -1) and (self.__restrictions[graph_tuple].find("u_turn") == -1):
                            include_num += 1
                            pass
                                
                        elif self.__restrictions[graph_tuple].find("only_") != -1:
                            include_num += 1
                            if gval["v"] == nval["u"]:
                                self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                            else:
                                pass
                    else:
                        if gval["v"] == nval["u"]:
                            self.__nxdual_rturns.append(append_to_dict(gkey, gval, nkey, nval))
                        

    @staticmethod
    def make_weakly_connected(graph_nx: nx.DiGraph) -> None:
        larges_ccomp = max(nx.strongly_connected_components(graph_nx), key=len)
        largest_graph_component = nx.DiGraph()
        largest_graph_component.add_nodes_from((n, graph_nx.nodes[n]) for n in larges_ccomp)
        for abc in graph_nx.edges.data():
            if (abc[0] in larges_ccomp) and (abc[1] in larges_ccomp):
                largest_graph_component.add_edge(abc[0], abc[1], len_m = abc[2]["len_m"])
        return largest_graph_component

    def make_graph_noturns(self) -> nx.DiGraph:
        self.__make_dual_noturns()
        nxgraph = nx.DiGraph()
        for val in self.__nxdual_noturns:
            nxgraph.add_edge(val["u"], val["v"], len_m = val["meters"])
        return nxgraph
    
    def make_graph_rturns(self) -> nx.DiGraph:
        self.__make_dual_withrturns()
        nxgraph = nx.DiGraph()
        for val in self.__nxdual_rturns:
            nxgraph.add_edge(val["u"], val["v"], len_m = val["meters"])
        return nxgraph

    def return_structure(self, restriicted_turns: bool) -> dict:
        if restriicted_turns:
            return self.__reverse_data_rturns
        else:
            return self.__reverse_data_noturns
        
    
def calculate_centralities(dual_graph_noturns: nx.DiGraph, dual_graph_turns: nx.DiGraph, graph_represent: nxGraph, shape_folder: str) -> None:
    bc_noturn = nx.betweenness_centrality(dual_graph_noturns, weight="meter_dist")
    cc_noturn = nx.closeness_centrality(dual_graph_noturns, distance="meter_dist")
    redir_matrix_noturn = graph_represent.return_structure(False)
    bc_list = []
    osmid_list = []
    highway_list = []
    oneway_list = []
    meter_dist_list = []
    geometry_list = []
    cc_list = []

    for key, val in redir_matrix_noturn.items():
        if key in bc_noturn:
            bc_list.append(bc_noturn[key])
            cc_list.append(cc_noturn[key])
            osmid_list.append(val["osm_id"])
            highway_list.append(val["highway"])
            oneway_list.append(val["oneway"])
            meter_dist_list.append(val["meter_dist"])
            geometry_list.append(val["geometry"])
            
    primal_gdf_noturns = gpd.GeoDataFrame({"osm_id": osmid_list, "bc": bc_list, "cc": cc_list, "highway": highway_list,
                                        "oneway": oneway_list, "meter_dist": meter_dist_list, "geometry": geometry_list})
    primal_gdf_noturns.set_crs(epsg=4326, inplace=True)
    primal_gdf_noturns.to_file(f"{shape_folder}/no_turns.shp")

    bc_list.clear()
    cc_list.clear()
    osmid_list.clear()
    highway_list.clear()
    oneway_list.clear()
    meter_dist_list.clear()
    geometry_list.clear()

    bc_turns = nx.betweenness_centrality(dual_graph_turns, weight="meter_dist")
    cc_turns = nx.closeness_centrality(dual_graph_turns, distance="meter_dist")
    redir_matrix = graph_represent.return_structure(True)

    for key, val in redir_matrix.items():
        if key in bc_turns:
            bc_list.append(bc_turns[key])
            cc_list.append(cc_turns[key])
            osmid_list.append(val["osm_id"])
            highway_list.append(val["highway"])
            oneway_list.append(val["oneway"])
            meter_dist_list.append(val["meter_dist"])
            geometry_list.append(val["geometry"])
            
    primal_gdf_turns = gpd.GeoDataFrame({"osm_id": osmid_list, "bc": bc_list, "cc": cc_list, "highway": highway_list,
                                        "oneway": oneway_list, "meter_dist": meter_dist_list, "geometry": geometry_list})
    primal_gdf_turns.set_crs(epsg=4326, inplace=True)
    primal_gdf_turns.to_file(f"{shape_folder}/with_turns.shp")

    difference_gdf = deepcopy(primal_gdf_noturns)
    diff_bc = []
    diff_cc = []

    for index, row in difference_gdf.iterrows():
        temp_layer = primal_gdf_turns.loc[primal_gdf_turns["geometry"] == row["geometry"]]
        try:
            diff_bc.append(row["bc"] - temp_layer["bc"].tolist()[0])
            diff_cc.append(row["cc"] - temp_layer["cc"].tolist()[0])
        except IndexError:
            diff_bc.append(0)
            diff_cc.append(0)

    difference_gdf["Diff_BC"] = diff_bc
    difference_gdf["Diff_CC"] = diff_cc

    difference_gdf.set_crs(epsg=4326, inplace=True)
    difference_gdf.to_file(f"{shape_folder}/difference.shp")

def zone_creation(nxgraph: nx.DiGraph, redir_matrix: dict, zone_size: list, osm_id: int, alpha, save_path: str) -> None:
    for k, v in redir_matrix.items():
        if v["osm_id"] == osm_id:
            zones = []
            for dist in zone_size:
                shortest_paths_lengths = nx.single_source_dijkstra_path_length(nxgraph, k, weight='len_m')
                accessible_nodes = [node for node, length in shortest_paths_lengths.items() if length <= dist]
                pseudo_node_geom = []
                for kv, vl in redir_matrix.items():
                    if kv in accessible_nodes:
                        for ptcoord in list(vl["geometry"].coords):
                            shpoint: shapely.Point = shapely.Point(ptcoord)
                            pseudo_node_geom.append(shpoint)
                try:
                    conc = shapely.concave_hull(shapely.unary_union(pseudo_node_geom), alpha)
                    print(conc)
                    zones.append(conc)
                except:
                    pass
            zones.reverse()
            print(zones)
            gdf = gpd.GeoDataFrame({"geometry": zones})
            gdf.set_crs(4326, inplace=True)
            gdf.to_file(save_path)