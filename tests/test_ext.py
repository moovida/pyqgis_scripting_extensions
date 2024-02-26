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

    def test_child_geometries(self):
        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly = HPolygon.fromCoords(coords1)
        self.assertEqual(poly.child_count(), 1)

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly1 = HPolygon.fromCoords(coords1)
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]
        poly2 = HPolygon.fromCoords(coords2)
        mp = HMultiPolygon([poly1, poly2])
        self.assertEqual(mp.child_count(), 2)
        child1 = mp.child(0)
        self.assertIsInstance(child1, HPolygon)

        ls = HLineString.fromCoords([(1, 2), (3, 4), (5, 6), (7, 8)])
        self.assertEqual(ls.child_count(), 1)

        ls1 = HLineString.fromCoords([(1, 2), (3, 4)])
        ls2 = HLineString.fromCoords([(5, 6), (7, 8)])
        mls = HMultiLineString([ls1, ls2])
        self.assertEqual(mls.child_count(), 2)

        p = HPoint(1, 2)
        self.assertEqual(p.child_count(), 1)

        p1 = HPoint(1, 2)
        p2 = HPoint(3, 4)
        mp = HMultiPoint([p1, p2])
        self.assertEqual(mp.child_count(), 2)

    def test_geometry_coordinates(self):
        coords = [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (7.0, 8.0), (1.0, 2.0)]
        poly = HPolygon.fromCoords(coords)
        self.assertEqual(poly.coordinates(), coords)

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8), (1, 2)]
        poly1 = HPolygon.fromCoords(coords1)
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16), (9, 10)]
        poly2 = HPolygon.fromCoords(coords2)
        mp = HMultiPolygon([poly1, poly2])
        self.assertEqual(mp.coordinates(), [coords1, coords2])

        coords = [(1, 2), (3, 4)]
        ls = HLineString.fromCoords(coords)
        self.assertEqual(ls.coordinates(), coords)

        coords1 = [(1, 2), (3, 4), (5, 6), (7, 8)]
        coords2 = [(9, 10), (11, 12), (13, 14), (15, 16)]
        mls = HMultiLineString.fromCoords([coords1, coords2])
        self.assertEqual(mls.coordinates(), [coords1, coords2])

        coords = [(1, 2), (3, 4)]
        mp = HMultiPoint.fromCoords(coords)
        self.assertEqual(mp.coordinates(), coords)

        p = HPoint(1, 2)
        self.assertEqual(p.coordinates(), [(1, 2)])


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

        # test difference
        diffGeom = self.g1.difference(self.g6)
        self.assertEquals(diffGeom.asWkt(), "Polygon ((0 5, 3 5, 3 3, 5 3, 5 0, 0 0, 0 5))")

    def test_buffer(self):
        bufferGeom = self.g1.buffer(1, 1, JOINSTYLE_BEVEL, ENDCAPSTYLE_SQUARE)
        print(bufferGeom.asWkt())
        self.assertEquals(bufferGeom.asWkt(), "Polygon ((0 -1, -1 0, -1 5, 0 6, 5 6, 6 5, 6 0, 5 -1, 0 -1))")


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




