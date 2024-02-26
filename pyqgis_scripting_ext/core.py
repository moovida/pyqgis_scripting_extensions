from qgis.core import *
from qgis.gui import *
from PyQt5.QtGui import QColor

ENDCAPSTYLE_ROUND = Qgis.EndCapStyle.Round
ENDCAPSTYLE_FLAT = Qgis.EndCapStyle.Flat
ENDCAPSTYLE_SQUARE = Qgis.EndCapStyle.Square

JOINSTYLE_ROUND = Qgis.JoinStyle.Round
JOINSTYLE_BEVEL = Qgis.JoinStyle.Bevel
JOINSTYLE_MITER = Qgis.JoinStyle.Miter


class HGeometry:
    def __init__(self, geometry: any):
        self.geometry = geometry
        self.qgsGeometry = QgsGeometry(self.geometry.clone())

    def __str__(self):
        return self.asWkt()
    
    def asWkt(self) -> str:
        return self.qgsGeometry.asWkt()
    
    def coordinates(self) -> list[tuple[float, float]]:
        if isinstance(self.geometry, HPoint):
            return [(self.x, self.y)]  
        elif isinstance(self.geometry, HLineString):
            return [(p.x, p.y) for p in self.geometry.points]
        elif isinstance(self.geometry, HPolygon):
            return [p.coordinates() for p in self.geometry.exterior_ring().points]
        elif isinstance(self.geometry, QgsPoint):
            return [(self.geometry.x(), self.geometry.y())]
        elif isinstance(self.geometry, QgsMultiPoint):
            return [(p.x(), p.y()) for p in self.geometry.vertices()]
        elif isinstance(self.geometry, QgsLineString):
            return [(p.x(), p.y()) for p in self.geometry.points()]
        elif isinstance(self.geometry, QgsMultiLineString):
            return [[(p.x(), p.y()) for p in line.points()] for line in self.geometry]
        elif isinstance(self.geometry, QgsPolygon):
            return [(p.x(), p.y()) for p in self.geometry.exteriorRing().points()]
        elif isinstance(self.geometry, QgsMultiPolygon):
            return [[(p.x(), p.y()) for p in ring.vertices()] for ring in self.geometry]
        else:
            raise ValueError(f"Unsupported geometry type: ${type(self.geometry)}")
        
    def area(self) -> float:
        return self.geometry.area()
    
    def length(self) -> float:
        return self.qgsGeometry.length()
    
    def centroid(self):
        return HGeometry.from_specialized(self.geometry.centroid())
    
    def distance(self, other) -> float:
        return self.qgsGeometry.distance(QgsGeometry(other.geometry.clone()))
    
    def bbox(self) -> list[float]:
        rect = self.geometry.boundingBox()
        return [rect.xMinimum(), rect.yMinimum(), rect.xMaximum(), rect.yMaximum()]
    
    def child_count(self) -> int:
        if isinstance(self.geometry, QgsPoint) or isinstance(self.geometry, QgsCurve) or isinstance(self.geometry, QgsLineString) or isinstance(self.geometry, QgsPolygon):
            return 1
        return self.geometry.childCount()
    
    def child(self, index: int):
        if isinstance(self.geometry, QgsPoint) or isinstance(self.geometry, QgsCurve) or isinstance(self.geometry, QgsLineString) or isinstance(self.geometry, QgsPolygon):
            return HGeometry.from_specialized(self.geometry)                
        childGeom = self.geometry.childGeometry(index)
        return HGeometry.from_specialized(childGeom)
    
    def intersects(self, other) -> bool:
        return self.qgsGeometry.intersects(QgsGeometry(other.geometry.clone()))
    
    def contains(self, other) -> bool:
        return self.qgsGeometry.contains(QgsGeometry(other.geometry.clone()))
    
    def touches(self, other) -> bool:
        return self.qgsGeometry.touches(QgsGeometry(other.geometry.clone()))
    
    def intersection(self, other):
        intersectionGeom =  self.qgsGeometry.intersection(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(intersectionGeom)
    
    def symdifference(self, other):
        symdifferenceGeom =  self.qgsGeometry.symDifference(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(symdifferenceGeom)
    
    def union(self, other):
        unionGeom =  self.qgsGeometry.combine(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(unionGeom)
    
    def difference(self, other):
        differenceGeom =  self.qgsGeometry.difference(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(differenceGeom)
    
    def buffer(self, distance: float, segments: int = 8, joinstyle: int = JOINSTYLE_ROUND, capstyle: int = ENDCAPSTYLE_ROUND):
        bufferGeom = self.qgsGeometry.buffer(distance, segments, capstyle, joinstyle, 1)
        return HGeometry.from_specialized(bufferGeom)
    
    

    @staticmethod
    def fromWkt(wkt: str):
        geom = QgsGeometry.fromWkt(wkt)
        return HGeometry.from_specialized(geom)

    @staticmethod
    def from_specialized(geom: any):
        if isinstance(geom, QgsGeometry):
            geom = geom.constGet()

        if not geom:
            return None
        if isinstance(geom, QgsPoint):
            return HPoint(qgs_point_geom=geom)
        elif isinstance(geom, QgsMultiPoint):
            return HMultiPoint(qgs_multipoint_geom=geom)
        elif isinstance(geom, QgsLineString):
            return HLineString(qgs_linestring_geom=geom)
        elif isinstance(geom, QgsMultiLineString):
            return HMultiLineString(qgs_multilinestring_geom=geom)
        elif isinstance(geom, QgsPolygon):
            poly = HPolygon(qgs_polygon_geom=geom)
            return poly
        elif isinstance(geom, QgsMultiPolygon):
            return HMultiPolygon(qgs_multipolygon_geom=geom)
        elif isinstance(geom, QgsCurve):
            return HLineString([HPoint(p.x(), p.y()) for p in geom.vertices()])
        else:
            raise ValueError(f"Unsupported geometry type: ${type(geom)}")
        
    

class HPoint(HGeometry):
    def __init__(self, x: float = None, y: float = None, qgs_point_geom: QgsPoint = None):
        if qgs_point_geom:
            self.x = qgs_point_geom.x()
            self.y = qgs_point_geom.y()
            self.geometry = qgs_point_geom
        else:
            self.x = x
            self.y = y
            self.geometry = QgsPoint(x, y)
        super().__init__(self.geometry)
    
    @classmethod
    def fromQgsPoint(cls, point: QgsPoint):
        return cls(point.x(), point.y())
    
    

class HMultiPoint(HGeometry):
    def __init__(self, points: list[HPoint] = None, qgs_multipoint_geom: QgsMultiPoint = None):
        if qgs_multipoint_geom:
            self.geometry = qgs_multipoint_geom
        else:
            self.geometry = QgsMultiPoint([QgsPoint(p.x, p.y) for p in points])
        super().__init__(self.geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[tuple[float, float]]):
        points = [HPoint(x, y) for x, y in coordinates]
        return cls(points = points)
    
class HLineString(HGeometry):
    def __init__(self, points: list[HPoint] = None, qgs_linestring_geom: QgsLineString = None):
        if qgs_linestring_geom:
            self.geometry = qgs_linestring_geom
        else:
            self.geometry = QgsLineString([QgsPoint(p.x, p.y) for p in points])
        super().__init__(self.geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[tuple[float, float]]):
        points = [HPoint(x, y) for x, y in coordinates]
        return cls(points)

class HMultiLineString(HGeometry):
    def __init__(self, linestrings: list[HLineString] = None, qgs_multilinestring_geom: QgsMultiLineString = None):
        if qgs_multilinestring_geom:
            self.geometry = qgs_multilinestring_geom
        else:
            self.linestrings = linestrings
            self.geometry = QgsMultiLineString()
            linesToAdd = [QgsLineString([QgsPoint(c[0], c[1]) for c in ls.coordinates()]) for ls in linestrings]
            for line in linesToAdd:
                self.geometry.addGeometry(line)
        super().__init__(self.geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[list[tuple[float, float]]]):
        linestrings = [HLineString.fromCoords(coords) for coords in coordinates]
        return cls(linestrings)

class HPolygon(HGeometry):
    def __init__(self, exterior_ring: HLineString = None, qgs_polygon_geom: QgsPolygon = None):
        if qgs_polygon_geom:
            self.geometry = qgs_polygon_geom
        else:
            coords = exterior_ring.coordinates()
            if coords[0] != coords[-1]:
                raise ValueError("The first and last point of the exterior ring must be equal")
            self.geometry = QgsPolygon(QgsLineString([QgsPoint(x, y) for x, y in coords]))
        super().__init__(self.geometry)

    @classmethod
    def fromCoords(cls, coords: list[tuple[float, float]]):
        exterior_ring = HLineString.fromCoords(coords)
        return cls(exterior_ring)

    def add_interior_ring(self, ring: HLineString):
        # check if ring has first and last point equal
        coords = ring.coordinates()
        if coords[0] != coords[-1]:
            raise ValueError("The first and last point of an interior ring must be equal")
        self.geometry.addInteriorRing(ring.geometry)

    def exterior_ring(self) -> HLineString:
        return HGeometry.from_specialized(self.geometry.exteriorRing())
    
    def interior_rings_count(self):
        return self.geometry.numInteriorRings()
    
    def interior_ring(self, index: int):
        return HGeometry.from_specialized(self.geometry.interiorRing(index))

class HMultiPolygon(HGeometry):
    def __init__(self, polygons: list[HPolygon] = None, qgs_multipolygon_geom: QgsMultiPolygon = None):
        if qgs_multipolygon_geom:
            self.geometry = qgs_multipolygon_geom
        else:
            qgsPolygons = [poly.geometry for poly in polygons]
            self.geometry = QgsMultiPolygon()
            for poly in qgsPolygons:
                self.geometry.addGeometry(poly)
        super().__init__(self.geometry)

    @classmethod
    def fromCoords(cls, coords: list[list[tuple[float, float]]]):
        polygons = [HPolygon.fromCoords(c) for c in coords]
        return cls(polygons)


class HMapCanvas():
    def __init__(self, iface: QgisInterface):
        if iface:
            self.iface = iface
            self.canvas = iface.mapCanvas()
        else:
            self.canvas = QgsMapCanvas()

    @staticmethod
    def new():
        return HMapCanvas(None)
    
    def set_extent(self, eswn: list[float]):
        self.canvas.setExtent(QgsRectangle(*eswn))

    def set_extent(self, extent: QgsRectangle):
        self.canvas.setExtent(extent)

    def add_geometry(self, geometry: HGeometry, color: str = 'red', width: float = 2.0):
        if isinstance(geometry, HPoint):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Point)
            r.addGeometry(QgsGeometry(geometry.geometry))
        elif isinstance(geometry, HLineString):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Line)
            r.setColor(QColor(color))
            r.setWidth(width)
            r.addGeometry(QgsGeometry(geometry.geometry))
        elif isinstance(geometry, HPolygon):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Polygon)
            r.addGeometry(QgsGeometry(geometry.geometry))
            fillColor = QColor(color)
            fillColor.setAlphaF(0.5)
            r.setColor(fillColor)
            strokeColor = QColor(color)
            r.setStrokeColor(strokeColor)
            r.setWidth(width)

    def remove_all_geometries(self):
        if self.iface:
            for item in self.iface.mapCanvas().scene().items():
                if isinstance(item, QgsRubberBand):
                    self.iface.mapCanvas().scene().removeItem(item)
        else:
            raise ValueError("Cannot remove geometries from a non-iface map canvas")
            

    def show(self):
        self.canvas.show()