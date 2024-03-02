import sys
sys.path.append("./")
import unittest
from pyqgis_scripting_ext.core import *



class TestPyQgisExt(unittest.TestCase):

    # before the tests create test geoms
    def setUp(self):
        self.g1 = HPolygon.fromCoords([[0, 0], [0, 5], [5, 5], [5, 0], [0, 0]])
        self.g2 = HPolygon.fromCoords([[5, 0], [5, 2], [7, 2], [7, 0], [5, 0]])
        self.g3 = HPoint(4, 1)
        self.g4 = HPoint(5, 4)
        self.g5 = HLineString.fromCoords([[1, 0], [1, 6]])
        self.g6 = HPolygon.fromCoords([[3, 3], [3, 6], [6, 6], [6, 3], [3, 3]])

    def test_HPoint(self):
        p = HPoint(1, 2)
        self.assertIsInstance(p, HPoint)
        self.assertIsInstance(p.geometry, QgsPoint)
        self.assertEqual(p.x, 1)
        self.assertEqual(p.y, 2)
        self.assertEqual(p.asWkt(), "Point (1 2)")

    def test_HLineString(self):
        p1 = HPoint(1, 2)
        p2 = HPoint(3, 4)
        ls = HLineString([p1, p2])
        self.assertIsInstance(ls, HLineString)
        self.assertIsInstance(ls.geometry, QgsLineString)
        self.assertEqual(ls.asWkt(), "LineString (1 2, 3 4)")
    
    def test_hlinestring_from_coordinates_list(self):
        ls = HLineString.fromCoords([(1, 2), (3, 4)])
        self.assertIsInstance(ls, HLineString)
        self.assertIsInstance(ls.geometry, QgsLineString)
        self.assertEqual(ls.asWkt(), "LineString (1 2, 3 4)")
        
        ls = HLineString.fromCoords([[1, 2], [3, 4]])
        self.assertIsInstance(ls, HLineString)
        self.assertIsInstance(ls.geometry, QgsLineString)
        self.assertEqual(ls.asWkt(), "LineString (1 2, 3 4)")

    def test_HPolygon(self):
        # coordinates list one
        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        ls1 = HLineString.fromCoords(coords1)

        poly = HPolygon(ls1)
        self.assertIsInstance(poly, HPolygon)
        self.assertIsInstance(poly.geometry, QgsPolygon)
        self.assertEqual(poly.asWkt(), "Polygon ((1 2, 3 4, 5 6, 7 8, 1 2))")

    def test_hpolygon_from_coords_with_interior_ring(self):
        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly = HPolygon.fromCoords(coords1)
        holeCoords = [(9, 10), (11, 12), (13, 14), (9, 10)]
        hole = HLineString.fromCoords(holeCoords)
        poly.add_interior_ring(hole)

        self.assertIsInstance(poly, HPolygon)
        self.assertIsInstance(poly.geometry, QgsPolygon)
        self.assertEqual(poly.interior_rings_count(), 1)
        self.assertEqual(poly.asWkt(), "Polygon ((1 2, 3 4, 5 6, 7 8, 1 2),(9 10, 11 12, 13 14, 9 10))")

    def test_HMultiPoint(self):
        p1 = HPoint(1, 2)
        p2 = HPoint(3, 4)
        mp = HMultiPoint([p1, p2])
        self.assertIsInstance(mp, HMultiPoint)
        self.assertIsInstance(mp.geometry, QgsMultiPoint)
        self.assertEqual(mp.asWkt(), "MultiPoint ((1 2),(3 4))")

    def test_HMultiPoint_from_coords(self):
        mp = HMultiPoint.fromCoords([(1, 2), (3, 4)])
        self.assertIsInstance(mp, HMultiPoint)
        self.assertIsInstance(mp.geometry, QgsMultiPoint)
        self.assertEqual(mp.asWkt(), "MultiPoint ((1 2),(3 4))")

    def test_HMultiLineString(self):
        ls1 = HLineString.fromCoords([(1, 2), (3, 4)])
        ls2 = HLineString.fromCoords([(5, 6), (7, 8)])
        mls = HMultiLineString([ls1, ls2])
        self.assertIsInstance(mls, HMultiLineString)
        self.assertIsInstance(mls.geometry, QgsMultiLineString)
        self.assertEqual(mls.asWkt(), "MultiLineString ((1 2, 3 4),(5 6, 7 8))")

    def test_HMultiLineString_from_coords(self):
        mls = HMultiLineString.fromCoords([[(1, 2), (3, 4)], [(5, 6), (7, 8)]])
        self.assertIsInstance(mls, HMultiLineString)
        self.assertIsInstance(mls.geometry, QgsMultiLineString)
        self.assertEqual(mls.asWkt(), "MultiLineString ((1 2, 3 4),(5 6, 7 8))")

    def test_HMultiPolygon(self):
        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly1 = HPolygon.fromCoords(coords1)
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]
        poly2 = HPolygon.fromCoords(coords2)
        mp = HMultiPolygon([poly1, poly2])
        self.assertIsInstance(mp, HMultiPolygon)
        self.assertIsInstance(mp.geometry, QgsMultiPolygon)
        self.assertEqual(mp.asWkt(), "MultiPolygon (((1 2, 3 4, 5 6, 7 8, 1 2)),((9 10, 11 12, 13 14, 15 16, 9 10)))")

    def test_HMultiPolygon_from_coords(self):
        mp = HMultiPolygon.fromCoords([[(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)], [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]])
        self.assertIsInstance(mp, HMultiPolygon)
        self.assertIsInstance(mp.geometry, QgsMultiPolygon)
        self.assertEqual(mp.asWkt(), "MultiPolygon (((1 2, 3 4, 5 6, 7 8, 1 2)),((9 10, 11 12, 13 14, 15 16, 9 10)))")

    def test_HGeometryCollection(self):
        p1 = HPoint(1, 2)
        p2 = HPoint(3, 4)
        mp = HMultiPoint([p1, p2])
        ls1 = HLineString.fromCoords([(1, 2), (3, 4)])
        ls2 = HLineString.fromCoords([(5, 6), (7, 8)])
        mls = HMultiLineString([ls1, ls2])
        gc = HGeometryCollection([mp, mls])
        self.assertIsInstance(gc, HGeometryCollection)
        self.assertIsInstance(gc.geometry, QgsGeometryCollection)
        self.assertEqual(gc.asWkt(), "GeometryCollection (MultiPoint ((1 2),(3 4)),MultiLineString ((1 2, 3 4),(5 6, 7 8)))")

    def test_child_geometries(self):
        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly = HPolygon.fromCoords(coords1)
        self.assertEqual(len(poly.geometries()), 1)

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly1 = HPolygon.fromCoords(coords1)
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]
        poly2 = HPolygon.fromCoords(coords2)
        mp = HMultiPolygon([poly1, poly2])

        subGeometries = mp.geometries()
        self.assertEqual(len(subGeometries), 2)
        child1 = subGeometries[0]
        self.assertIsInstance(child1, HPolygon)

        ls = HLineString.fromCoords([(1, 2), (3, 4), (5, 6), (7, 8)])
        self.assertEqual(len(ls.geometries()), 1)

        ls1 = HLineString.fromCoords([(1, 2), (3, 4)])
        ls2 = HLineString.fromCoords([(5, 6), (7, 8)])
        mls = HMultiLineString([ls1, ls2])
        self.assertEqual(len(mls.geometries()), 2)

        p = HPoint(1, 2)
        self.assertEqual(len(p.geometries()), 1)

        p1 = HPoint(1, 2)
        p2 = HPoint(3, 4)
        mp = HMultiPoint([p1, p2])
        self.assertEqual(len(mp.geometries()), 2)

    def test_geometry_coordinates(self):
        coords = [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (7.0, 8.0), (1.0, 2.0)]
        poly = HPolygon.fromCoords(coords)
        self.assertTrue(self.coordsEqual(poly.coordinates(), coords))

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly1 = HPolygon.fromCoords(coords1)
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]
        poly2 = HPolygon.fromCoords(coords2)
        mp = HMultiPolygon([poly1, poly2])
        self.assertTrue(self.coordsEqual2(mp.coordinates(), [coords1, coords2]))

        coords = [(1, 2), (3, 4)]
        ls = HLineString.fromCoords(coords)
        self.assertTrue(self.coordsEqual(ls.coordinates(), coords))

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8)]
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16)]
        mls = HMultiLineString.fromCoords([coords1, coords2])
        self.assertTrue(self.coordsEqual2(mls.coordinates(), [coords1, coords2]))

        coords = [(1, 2), (3, 4)]
        mp = HMultiPoint.fromCoords(coords)
        self.assertTrue(self.coordsEqual(mp.coordinates(), coords))

        p = HPoint(1, 2)
        self.assertTrue(self.coordsEqual(p.coordinates(), [(1, 2)]))

    def coordsEqual(self, coords1:list[float], coords2:list[float]):
        for i in range(len(coords1)):
            if coords1[i][0] != coords2[i][0] or coords1[i][1] != coords2[i][1]:
                return False
        return True

    def coordsEqual2(self, coords1:list[list[float]], coords2:list[list[float]]):
        for i in range(len(coords1)):
            for j in range(len(coords1[i])):
                if coords1[i][j][0] != coords2[i][j][0] or coords1[i][j][1] != coords2[i][j][1]:
                    return False
        return True


    def test_geometry_from_wkt(self):
        wkt = "POINT (1 2)"
        p = HGeometry.fromWkt(wkt)
        self.assertIsInstance(p, HPoint)

        wkt = "MULTIPOINT ((1 2),(3 4))"
        mp = HGeometry.fromWkt(wkt)
        self.assertIsInstance(mp, HMultiPoint)

        wkt = "LINESTRING (1 2, 3 4)"
        ls = HGeometry.fromWkt(wkt)
        self.assertIsInstance(ls, HLineString)

        wkt = "MULTILINESTRING ((1 2, 3 4),(5 6, 7 8))"
        mls = HGeometry.fromWkt(wkt)
        self.assertIsInstance(mls, HMultiLineString)

        wkt = "POLYGON ((1 2, 3 4, 5 6, 7 8, 1 2))"
        poly = HGeometry.fromWkt(wkt)
        self.assertIsInstance(poly, HPolygon)

        wkt = "MULTIPOLYGON (((1 2, 3 4, 5 6, 7 8, 1 2)),((9 10, 11 12, 13 14, 15 16, 9 10)))"
        mp = HGeometry.fromWkt(wkt)
        self.assertIsInstance(mp, HMultiPolygon)

    def test_area_and_length(self):
        self.assertEqual(self.g1.area(), 25)
        self.assertEqual(self.g1.length(), 20)

        self.assertEqual(self.g2.area(), 4)
        self.assertEqual(self.g2.length(), 8)

        self.assertEqual(self.g3.area(), 0)
        self.assertEqual(self.g3.length(), 0)

        self.assertEqual(self.g4.area(), 0)
        self.assertEqual(self.g4.length(), 0)

        self.assertEqual(self.g5.area(), 0)
        self.assertEqual(self.g5.length(), 6)
        
        self.assertEqual(self.g6.area(), 9)
        self.assertEqual(self.g6.length(), 12)

    def test_boundingbox(self):
        self.assertEqual(self.g1.bbox(), [0, 0, 5, 5])
        self.assertEqual(self.g2.bbox(), [5, 0, 7, 2])
        self.assertEqual(self.g3.bbox(), [4, 1, 4, 1])
        self.assertEqual(self.g4.bbox(), [5, 4, 5, 4])
        self.assertEqual(self.g5.bbox(), [1, 0, 1, 6])
        self.assertEqual(self.g6.bbox(), [3, 3, 6, 6])

    def test_distance(self):
        self.assertEqual(self.g4.distance(self.g5), 4)
        self.assertEqual(self.g5.distance(self.g6), 2)

    def test_predicate_methods(self):
        self.assertTrue(self.g1.intersects(self.g2))
        self.assertTrue(self.g1.intersects(self.g3))
        self.assertTrue(self.g1.intersects(self.g4))
        self.assertTrue(self.g1.intersects(self.g5))
        self.assertTrue(self.g1.intersects(self.g6))

        self.assertTrue(self.g1.touches(self.g2))
        self.assertFalse(self.g1.touches(self.g3))
        self.assertTrue(self.g1.touches(self.g4))
        self.assertFalse(self.g1.touches(self.g5))
        self.assertFalse(self.g1.touches(self.g6))

        self.assertFalse(self.g1.contains(self.g2))
        self.assertTrue(self.g1.contains(self.g3))
        self.assertFalse(self.g1.contains(self.g4))
        self.assertFalse(self.g1.contains(self.g5))
        self.assertFalse(self.g1.contains(self.g6))

    def test_functions(self):
        # test intersection
        intersectionGeom = self.g1.intersection(self.g6)
        self.assertEquals(intersectionGeom.asWkt(), "Polygon ((5 5, 5 3, 3 3, 3 5, 5 5))")

        # test symdiff
        symdiffGeom = self.g1.symdifference(self.g6)
        self.assertEquals(symdiffGeom.asWkt(), "MultiPolygon (((0 5, 3 5, 3 3, 5 3, 5 0, 0 0, 0 5)),((5 5, 3 5, 3 6, 6 6, 6 3, 5 3, 5 5)))")

        # test union
        unionGeom = self.g1.union(self.g6)
        self.assertEquals(unionGeom.asWkt(), "Polygon ((0 5, 3 5, 3 6, 6 6, 6 3, 5 3, 5 0, 0 0, 0 5))")

        unionGeom = self.g1 + self.g6
        self.assertEquals(unionGeom.asWkt(), "Polygon ((0 5, 3 5, 3 6, 6 6, 6 3, 5 3, 5 0, 0 0, 0 5))")

        # test difference
        diffGeom = self.g1.difference(self.g6)
        self.assertEquals(diffGeom.asWkt(), "Polygon ((0 5, 3 5, 3 3, 5 3, 5 0, 0 0, 0 5))")

        diffGeom = self.g1 - self.g6
        self.assertEquals(diffGeom.asWkt(), "Polygon ((0 5, 3 5, 3 3, 5 3, 5 0, 0 0, 0 5))")

    def test_buffer(self):
        bufferGeom = self.g1.buffer(1, 1, JOINSTYLE_BEVEL, ENDCAPSTYLE_SQUARE)
        self.assertEquals(bufferGeom.asWkt(), "Polygon ((0 -1, -1 0, -1 5, 0 6, 5 6, 6 5, 6 0, 5 -1, 0 -1))")

    def test_transform(self):
        crsHelper = HCrs()
        crsHelper.from_srid(4326)
        crsHelper.to_srid(32632)

        point4326 = HPoint(11, 46)
        point32632 = crsHelper.transform(point4326)

        self.assertEquals(int(point32632.x), 654863)
        self.assertEquals(int(point32632.y), 5095992)

        point4326back = crsHelper.transform(point32632, inverse=True)
        self.assertEquals(round(point4326back.x), 11)
        self.assertEquals(round(point4326back.y), 46)

    def test_vectorlayer_read(self):
        # read the vector layer
        layer = HVectorLayer.open("tests/ne.gpkg", "ne_10m_populated_places")
        
        self.assertEquals(layer.prjcode, "EPSG:4326")
        self.assertEquals(layer.size(), 58)

        fields = layer.fields
        self.assertEquals(fields["SCALERANK"], "Integer64")
        self.assertEquals(fields["ADM0_A3"], "String")

        # get all features
        features = layer.features()
        self.assertEquals(len(features), 58)

        # get features with filter
        features = layer.features("POP_MAX > 3300000")
        self.assertEquals(len(features), 1)
        cityName = features[0].attributes[layer.field_index('NAME')]
        self.assertEquals(cityName, "Rome")

        # get features with bbox filter
        features = layer.features(bbox=[12, 41, 13, 42])
        self.assertEquals(len(features), 2)

        # get features in polygon geometry
        poly = HPolygon.fromCoords([[12, 41], [13, 41], [13, 42], [12, 42], [12, 41]])
        features = layer.features(geometryfilter=poly)
        self.assertEquals(len(features), 2)

        self.assertIsInstance(features[0].geometry, HPoint)

    def test_vectorlayer_sublayer(self):
        # read the vector layer
        layer = HVectorLayer.open("tests/ne.gpkg", "ne_10m_populated_places")
        minX = 10.5
        maxX = 11.5
        minY= 45.2
        maxY = 46.8
        intersectionPolygon = HPolygon.fromCoords([[minX, minY], [maxX, minY], [maxX, maxY], [minX, maxY], [minX, minY]])

        # complete attributes sublayer
        subLayer = layer.sub_layer(intersectionPolygon, "sublayer_test")
        self.assertEquals(subLayer.size(), 3)
        index = subLayer.field_index('SCALERANK')
        self.assertEquals(index, 1)
        index = subLayer.field_index('NAME')
        self.assertEquals(index, 5)

        # only extract the NAME attribute
        subLayer = layer.sub_layer(intersectionPolygon, "sublayer_test", ["NAME"])
        self.assertEquals(subLayer.size(), 3)
        index = subLayer.field_index('SCALERANK')
        self.assertEquals(index, -1)
        index = subLayer.field_index('NAME')
        self.assertEquals(index, 0)
        names = [f.attributes[index] for f in subLayer.features()]
        self.assertIn("Bolzano", names)
        self.assertIn("Trento", names)
        self.assertIn("Verona", names)


        # geopackagePath = "/storage/data/natural_earth_vector.gpkg"
        # countriesName = "ne_10m_admin_0_countries"
        # riversName = "ne_10m_rivers_lake_centerlines_scale_rank"
        # riversLayer = HVectorLayer.open(geopackagePath, riversName)
        # countriesLayer = HVectorLayer.open(geopackagePath, countriesName)
        
        # countriesLayer.subset_filter("NAME='Italy'")
        # italyFeatures = countriesLayer.features("NAME='Italy'")

        # riversItalyLayer = riversLayer.sub_layer(italyFeatures[0].geometry, "rivers_italy", ['scalerank', 'name'])
        # features = riversItalyLayer.features()
        # print("rivers", len(features))
        # for f in features:
        #     print(f)




    def test_vectorlayer_write(self):
        fields = {
            "id": "Integer",
            "name": "String",
        }
        layer = HVectorLayer.new("test", "Point", "EPSG:4326", fields)
        layer.add_feature(HPoint(-122.42, 37.78), [1, "San Francisco"])
        layer.add_feature(HPoint(-73.98, 40.47), [2, "New York"])

        self.assertEquals(len(layer.features()), 2)
        self.assertEquals(len(layer.features("name='New York'")), 1)
        self.assertEquals(layer.prjcode, "EPSG:4326")

        # test geopackage
        error = layer.dump_to_gpkg("tests/test.gpkg", overwrite=True)
        if(error):
            print(error )
        self.assertIsNone(error)
        self.assertTrue(os.path.exists("tests/test.gpkg"))
        os.remove("tests/test.gpkg")
        if os.path.exists("tests/test.gpkg-shm"):
            os.remove("tests/test.gpkg-shm")
        if os.path.exists("tests/test.gpkg-wal"):
            os.remove("tests/test.gpkg-wal")
        
        # test shp
        error = layer.dump_to_shp("tests/test.shp")
        if(error):
            print(error )
        self.assertIsNone(error)
        self.assertTrue(os.path.exists("tests/test.shp"))
        os.remove("tests/test.shp")
        if os.path.exists("tests/test.shx"):
            os.remove("tests/test.shx")
        if os.path.exists("tests/test.dbf"):
            os.remove("tests/test.dbf")
        if os.path.exists("tests/test.prj"):
            os.remove("tests/test.prj")
        if os.path.exists("tests/test.cpg"):
            os.remove("tests/test.cpg")


        


    def test_style_composing(self):
        pointStyle = HMarker("circle", 10) + HFill("green") + HStroke("black", 1)
        self.assertEquals(pointStyle.type, "point")
        pointStyle = HFill("green") + HStroke("black", 1) + HMarker("circle", 10)
        self.assertEquals(pointStyle.type, "point")
        pointStyle =  HStroke("black", 1) + HMarker("circle", 10) + HFill("green")
        self.assertEquals(pointStyle.type, "point")

        lineStyle = HStroke("blue", 2)
        self.assertEquals(lineStyle.type, "line")

        polyStyle = HFill("green") + HStroke("black", 1)
        self.assertEquals(polyStyle.type, "polygon")

        labelProperties = {
            "font": "Arial", 
            "color": 'black', 
            "size": 10.0, 
            "field": "dummy", 
            "xoffset": 10.0, 
            "yoffset": 10.0
        }

        point_with_label = HMarker("circle", 10) + HFill("green") + HStroke("black", 1) + HLabel(**labelProperties) + HHalo("white", 1)
        self.assertEquals(point_with_label.type, "point")
        self.assertEquals(point_with_label.properties['label_font'], "Arial")
        self.assertEquals(point_with_label.properties['halo_color'], "white")
        

        
        


    # def test_mapcanvas(self):
    #     canvas = HMapCanvas.new()

    #     # create a geometry per type
    #     p = HPoint(1, 2)
    #     ls = HLineString.fromCoords([(1, 2), (3, 4)])
    #     poly = HPolygon.fromCoords([(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)])
        
    #     # add the geometry to the map canvas with color and width
    #     canvas.add_geometry(p, "red", 2)
    #     canvas.add_geometry(ls, "blue", 2)
    #     canvas.add_geometry(poly, "green", 2)

    #     # set the extent of the map canvas
    #     canvas.set_extent([0, 0, 10, 10])

    #     # show
    #     canvas.show()




