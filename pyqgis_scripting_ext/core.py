from qgis.core import *
from qgis.gui import *
import os
from PyQt5.QtGui import QColor

ENDCAPSTYLE_ROUND = Qgis.EndCapStyle.Round
ENDCAPSTYLE_FLAT = Qgis.EndCapStyle.Flat
ENDCAPSTYLE_SQUARE = Qgis.EndCapStyle.Square

JOINSTYLE_ROUND = Qgis.JoinStyle.Round
JOINSTYLE_BEVEL = Qgis.JoinStyle.Bevel
JOINSTYLE_MITER = Qgis.JoinStyle.Miter


class HGeometry:
    def __init__(self, geometry: QgsAbstractGeometry):
        self.geometry = geometry
        self.qgs_geometry = None

    def get_qgs_geometry(self) -> QgsGeometry:
        if not self.qgs_geometry:
            self.qgs_geometry = QgsGeometry(self.geometry.clone())
        return self.qgs_geometry

    def __str__(self):
        return self.asWkt()
    
    def asWkt(self) -> str:
        return self.get_qgs_geometry().asWkt()
    
    def coordinates(self) -> list[tuple[float, float]]:
        """
        Get the coordinates of the geometry.
        """
        if isinstance(self.geometry, HPoint):
            return [(self.geometry.x(), self.geometry.y(), self.geometry.z())]  
        elif isinstance(self.geometry, HLineString):
            return [(p.x, p.y, p.z) for p in self.geometry.geometry.points()]
        elif isinstance(self.geometry, HPolygon):
            return [p.coordinates() for p in self.geometry.exterior_ring().geometry.points()]
        elif isinstance(self.geometry, QgsPoint):
            return [(self.geometry.x(), self.geometry.y(), self.geometry.z())]
        elif isinstance(self.geometry, QgsMultiPoint):
            return [(p.x(), p.y(), p.z()) for p in self.geometry.vertices()]
        elif isinstance(self.geometry, QgsLineString):
            return [(p.x(), p.y(), p.z()) for p in self.geometry.points()]
        elif isinstance(self.geometry, QgsMultiLineString):
            return [[(p.x(), p.y(), p.z()) for p in line.points()] for line in self.geometry]
        elif isinstance(self.geometry, QgsPolygon):
            return [(p.x(), p.y(), p.z()) for p in self.geometry.exteriorRing().points()]
        elif isinstance(self.geometry, QgsMultiPolygon):
            return [[(p.x(), p.y(), p.z()) for p in ring.exteriorRing().points()] for ring in self.geometry]
        elif isinstance(self.geometry, QgsGeometryCollection):
            return [HGeometry.from_specialized(g).coordinates() for g in self.geometry]
        else:
            raise ValueError(f"Unsupported geometry type: ${type(self.geometry)}")
        
    def area(self) -> float:
        return self.geometry.area()
    
    def length(self) -> float:
        return self.get_qgs_geometry().length()
    
    def centroid(self):
        return HGeometry.from_specialized(self.geometry.centroid())
    
    def distance(self, other) -> float:
        return self.get_qgs_geometry().distance(QgsGeometry(other.geometry.clone()))
    
    def bbox(self) -> list[float]:
        rect = self.geometry.boundingBox()
        return [rect.xMinimum(), rect.yMinimum(), rect.xMaximum(), rect.yMaximum()]
    
    def geometries(self):
        """
        Get the child geometries of the geometry.

        If the geometry is a point, line or polygon, it returns a list with the geometry itself.
        """
        if isinstance(self.geometry, QgsPoint) or isinstance(self.geometry, QgsCurve) or isinstance(self.geometry, QgsLineString) or isinstance(self.geometry, QgsPolygon):
            return [HGeometry.from_specialized(self.geometry)]
        count = self.geometry.childCount()
        return [HGeometry.from_specialized(self.geometry.childGeometry(i)) for i in range(count)]
    
    def intersects(self, other) -> bool:
        return self.get_qgs_geometry().intersects(QgsGeometry(other.geometry.clone()))
    
    def contains(self, other) -> bool:
        return self.get_qgs_geometry().contains(QgsGeometry(other.geometry.clone()))
    
    def touches(self, other) -> bool:
        return self.get_qgs_geometry().touches(QgsGeometry(other.geometry.clone()))
    
    def overlaps(self, other) -> bool:
        return self.get_qgs_geometry().overlaps(QgsGeometry(other.geometry.clone()))
    
    def crosses(self, other) -> bool:
        return self.get_qgs_geometry().crosses(QgsGeometry(other.geometry.clone()))
    
    def intersection(self, other):
        intersectionGeom =  self.get_qgs_geometry().intersection(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(intersectionGeom)
    
    def symdifference(self, other):
        symdifferenceGeom =  self.get_qgs_geometry().symDifference(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(symdifferenceGeom)
    
    def union(self, other):
        unionGeom =  self.get_qgs_geometry().combine(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(unionGeom)

    def __add__(self, other):
        return self.union(other)
    
    def difference(self, other):
        differenceGeom =  self.get_qgs_geometry().difference(QgsGeometry(other.geometry.clone()))
        return HGeometry.from_specialized(differenceGeom)
    
    def __sub__(self, other):
        return self.difference(other)
    
    def buffer(self, distance: float, segments: int = 8, joinstyle: int = JOINSTYLE_ROUND, capstyle: int = ENDCAPSTYLE_ROUND):
        bufferGeom = self.get_qgs_geometry().buffer(distance, segments, capstyle, joinstyle, 1)
        return HGeometry.from_specialized(bufferGeom)
    
    def centroid(self):
        return HGeometry.from_specialized(self.get_qgs_geometry().centroid())
    
    def convex_hull(self):
        return HGeometry.from_specialized(self.get_qgs_geometry().convexHull())
    
    @staticmethod
    def fromWkt(wkt: str):
        geom = QgsGeometry.fromWkt(wkt)
        return HGeometry.from_specialized(geom)

    @staticmethod
    def from_specialized(geom: any):
        """
        Create a HGeometry object from a specialized geometry object (Qgs*).
        """
        if isinstance(geom, QgsGeometry):
            geom = geom.constGet()
        geom = geom.clone()

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
            return HLineString([HPoint(p.x(), p.y(), p.z()) for p in geom.vertices()])
        elif isinstance(geom, QgsGeometryCollection):
            return HGeometryCollection(qgs_geometrycollection_geom=geom)
        else:
            raise ValueError(f"Unsupported geometry type: ${type(geom)}")
        
    

class HPoint(HGeometry):
    def __init__(self, x: float = None, y: float = None, z: float = None , qgs_point_geom: QgsPoint = None):
        if qgs_point_geom:
            self.x = qgs_point_geom.x()
            self.y = qgs_point_geom.y()
            self.z = qgs_point_geom.z()
            geometry = qgs_point_geom
        else:
            self.x = x
            self.y = y
            self.z = z
            geometry = QgsPoint(x, y, z)
        super().__init__(geometry)
        
    
    def x(self) -> float:
        return self.geometry.x()
    
    def y(self) -> float:
        return self.geometry.y()
    
    def z(self) -> float:
        return self.geometry.z()
    
    

class HMultiPoint(HGeometry):
    def __init__(self, points: list[HPoint] = None, qgs_multipoint_geom: QgsMultiPoint = None):
        if qgs_multipoint_geom:
            geometry = qgs_multipoint_geom
        else:
            geometry = QgsMultiPoint([QgsPoint(p.x, p.y, p.z) for p in points])
        super().__init__(geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[tuple[float, float]]):
        points = [HPoint(x, y) for x, y in coordinates]
        return cls(points = points)
    
class HLineString(HGeometry):
    def __init__(self, points: list[HPoint] = None, qgs_linestring_geom: QgsLineString = None):
        if qgs_linestring_geom:
            geometry = qgs_linestring_geom
        else:
            geometry = QgsLineString([QgsPoint(p.x, p.y) for p in points])
        super().__init__(geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[tuple[float, float]]):
        points = [HPoint(x, y) for x, y in coordinates]
        return cls(points)

class HMultiLineString(HGeometry):
    def __init__(self, linestrings: list[HLineString] = None, qgs_multilinestring_geom: QgsMultiLineString = None):
        if qgs_multilinestring_geom:
            geometry = qgs_multilinestring_geom
        else:
            geometry = QgsMultiLineString()
            linesToAdd = [QgsLineString([QgsPoint(c[0], c[1]) for c in ls.coordinates()]) for ls in linestrings]
            for line in linesToAdd:
                geometry.addGeometry(line)
        super().__init__(geometry)

    @classmethod
    def fromCoords(cls, coordinates: list[list[tuple[float, float]]]):
        linestrings = [HLineString.fromCoords(coords) for coords in coordinates]
        return cls(linestrings)

class HPolygon(HGeometry):
    def __init__(self, exterior_ring: HLineString = None, qgs_polygon_geom: QgsPolygon = None):
        if qgs_polygon_geom:
            geometry = qgs_polygon_geom
        else:
            coords = exterior_ring.coordinates()
            if coords[0][0] != coords[-1][0] or coords[0][1] != coords[-1][1]:
                raise ValueError("The first and last point of the exterior ring must be equal")
            geometry = QgsPolygon(QgsLineString([QgsPoint(x, y, z) for x, y, z in coords]))
        super().__init__(geometry)

    @classmethod
    def fromCoords(cls, coords: list[tuple[float, float]]):
        exterior_ring = HLineString.fromCoords(coords)
        return cls(exterior_ring)

    def add_interior_ring(self, ring: HLineString):
        # check if ring has first and last point equal
        coords = ring.coordinates()
        if coords[0][0] != coords[-1][0] or coords[0][1] != coords[-1][1]:
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
            geometry = qgs_multipolygon_geom
        else:
            qgsPolygons = [poly.geometry for poly in polygons]
            geometry = QgsMultiPolygon()
            for poly in qgsPolygons:
                geometry.addGeometry(poly)
        super().__init__(geometry)

    @classmethod
    def fromCoords(cls, coords: list[list[tuple[float, float]]]):
        polygons = [HPolygon.fromCoords(c) for c in coords]
        return cls(polygons)

class HGeometryCollection(HGeometry):
    def __init__(self, geometries: list[HGeometry] = None, qgs_geometrycollection_geom: QgsGeometryCollection = None):
        if qgs_geometrycollection_geom:
            geometry = qgs_geometrycollection_geom
        else:
            geometry = QgsGeometryCollection()
            for geom in geometries:
                geometry.addGeometry(geom.geometry)
        super().__init__(geometry)


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
        self.canvas.setExtent(QgsRectangle(eswn[0], eswn[1], eswn[2], eswn[3]))

    # def set_extent(self, extent: QgsRectangle):
    #     self.canvas.setExtent(extent)

    def add_geometry(self, geometry: HGeometry, color: str = 'red', width: float = 2.0):
        if isinstance(geometry, HPoint):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Point)
            r.addGeometry(QgsGeometry(geometry.geometry.clone()))
            r.setColor(QColor(color))
            r.setWidth(width)
        elif isinstance(geometry, HMultiPoint):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Point)
            for point in geometry.geometries():
                r.addGeometry(QgsGeometry(point.geometry.clone()))
            r.setColor(QColor(color))
            r.setWidth(width)
        elif isinstance(geometry, HLineString):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Line)
            r.setColor(QColor(color))
            r.setWidth(width)
            r.addGeometry(QgsGeometry(geometry.geometry.clone()))
        elif isinstance(geometry, HMultiLineString):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Line)
            r.setColor(QColor(color))
            r.setWidth(width)
            for line in geometry.geometries():
                r.addGeometry(QgsGeometry(line.geometry.clone()))
        elif isinstance(geometry, HPolygon):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Polygon)
            r.addGeometry(QgsGeometry(geometry.geometry.clone()))
            fillColor = QColor(color)
            fillColor.setAlphaF(0.5)
            r.setColor(fillColor)
            strokeColor = QColor(color)
            r.setStrokeColor(strokeColor)
            r.setWidth(width)
        elif isinstance(geometry, HMultiPolygon):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Polygon)
            for poly in geometry.geometries():
                r.addGeometry(QgsGeometry(poly.geometry.clone()))
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


class HCrs:
    def __init__(self):
        self.fromCrs = None
        self.toCrs = None
        self.crsTransf = None

    def from_srid(self, srid):
        if isinstance(srid, str):
            epsg = srid
        elif isinstance(srid, int):
            epsg = f"EPSG:{srid}"
        else:
            raise ValueError("The SRID must be a string or an integer")
        self.fromCrs = QgsCoordinateReferenceSystem(epsg)

    def to_srid(self, srid):
        if isinstance(srid, str):
            epsg = srid
        elif isinstance(srid, int):
            epsg = f"EPSG:{srid}"
        else:
            raise ValueError("The SRID must be a string or an integer")
        self.toCrs = QgsCoordinateReferenceSystem(epsg)

    def from_crs(self, crs: QgsCoordinateReferenceSystem):
        self.fromCrs = crs

    def to_crs(self, crs: QgsCoordinateReferenceSystem):
        self.toCrs = crs

    def transform(self, geometry: HGeometry, inverse:bool = False) -> HGeometry:
        """
        Transform a geometry from the fromCrs to the toCrs.
        """
        if not self.fromCrs or not self.toCrs:
            raise ValueError("The from and to CRS must be set")
        if not self.crsTransf:
            self.crsTransf = QgsCoordinateTransform(self.fromCrs, self.toCrs, QgsProject.instance())
        qgsGeometry = QgsGeometry(geometry.geometry.clone())
        if inverse:
            QgsGeometry.transform(qgsGeometry, self.crsTransf, direction=Qgis.TransformDirection.Reverse)
        else:
            QgsGeometry.transform(qgsGeometry, self.crsTransf)
        return HGeometry.from_specialized(qgsGeometry)
    
    def transfrom(self, x:float, y:float) -> tuple[float, float]:
        """
        Transform a point from the fromCrs to the toCrs.
        """
        if not self.fromCrs or not self.toCrs:
            raise ValueError("The from and to CRS must be set")
        if not self.crsTransf:
            self.crsTransf = QgsCoordinateTransform(self.from_crs, self.to_crs, QgsProject.instance())
        qgsPointXY = self.crsTransf.transform(x, y)
        return (qgsPointXY.x(), qgsPointXY.y())

    @staticmethod
    def current_crs():
        """
        Get the CRS of the current project.
        """
        return QgsProject.instance().crs()

class HFeature:
    def __init__(self, feature: QgsFeature):
        self.feature = feature
        self.geometry = HGeometry.from_specialized(self.feature.geometry())
        self.attributes = self.feature.attributeMap()

    @staticmethod
    def create(geometry: HGeometry, attributes: dict[str, any]):
        """
        Create a new feature based on the givengeometry and attributes.

        The attributes must be a dictionary with the field names as keys 
        and the same order as the schema they will be written to..
        """
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry(geometry.geometry.clone()))
        feature.setAttributes(attributes)
        return HFeature(feature)

class HStyle:
    def __init__(self, properties: dict, type:str = "line"):
        self.properties = properties
        if not type:
            type = "line"
        if not type in ["point", "line", "polygon"]:
            raise ValueError("The style type must be point, line or polygon")
        self.type = type

    def __add__(self, other):
        # new style with the properties of both styles
        properties = {}
        properties.update(self.properties)
        properties.update(other.properties)

        type = "line"
        if (self.type == "point" or other.type == "point"):
            type = "point"
        elif (self.type == "polygon" or other.type == "polygon"):
            type = "polygon"

        return HStyle(properties, type)
    
    def __str__(self) -> str:
        string = "type: " + self.type + "\n"
        for key, value in self.properties.items():
            string += f"\t{key}: {value}\n"
        return string

class HMarker(HStyle):
    def __init__(self, name: str = "square", size: float = 5.0, angle: float = 0.0):
        self.name = name
        self.size = size
        self.angle = angle
        properties = {
            "marker_name": name,
            "marker_size": size,
            "marker_angle": angle
        }
        super().__init__(properties)
        self.type = "point"

class HStroke(HStyle):
    def __init__(self, color: str = 'red', width: float = 2.0):
        self.color = color
        self.width = width
        properties = {
            "stroke_color": color,
            "stroke_width": width,
        }
        super().__init__(properties)

class HFill(HStyle):
    def __init__(self, color: str = 'red'):
        self.color = color
        properties = {
            "fill_color": color,
        }
        super().__init__(properties)
        if self.type == "line":
            self.type = "polygon"

class HLabel(HStyle):
    def __init__(self, font: str = "Arial", color: str = 'black', size: float = 10.0, field:str = None, \
                 xoffset: float = 0.0, yoffset: float = 0.0):    
        self.font = font
        self.color = color
        self.size = size
        self.field = field
        self.xoffset = xoffset
        self.yoffset = yoffset
        properties = {
            "label_font": font,
            "label_color": color,
            "label_size": size,
            "label_field": field,
            "label_xoffset": xoffset,
            "label_yoffset": yoffset,
        }
        super().__init__(properties)

class HVectorLayer:
    def __init__(self, layer: QgsVectorLayer, isReadOnly: bool):
        self.layer = layer
        self.isReadOnly = isReadOnly

    @staticmethod
    def open(path: str, table_name: str):
        """
        Open a vector layer from a file path.
        
        If the file is a geopackage, the table name must be specified.
        """
        # if it is a geopackage, the table name must be specified
        if path.endswith(".gpkg"):
            uri = f"{path}|layername={table_name}"
            layer = QgsVectorLayer(uri, table_name, "ogr")
        elif path.endswith(".shp"):
            # get file name without extension from path
            layer_name = os.path.splitext(os.path.basename(path))[0]
            layer = QgsVectorLayer(path, layer_name, "ogr")
        else:
            raise ValueError(f"Unsupported file format:{path}")
        return HVectorLayer(layer, True)
    
    @staticmethod
    def new(name: str, geometry_type: str, srid: int, fields: dict[str, str]):
        """
        Create a new vector layer in memory.
        """
        definition = f"{geometry_type}?crs=epsg:{srid}"
        for fname, type in fields.items():
            definition += f"&field={fname}:{type}"
        layer = QgsVectorLayer(definition, name, "memory")
        return HVectorLayer(layer, False)

    def srid(self) -> str:
        """
        Get the SRID of the layer.

        It is of the form "EPSG:xxxx"
        """
        return self.layer.crs().authid()
    
    def fields(self) -> dict:
        """
        Get the field namnes and types of the layer.
        """
        fieldsDict = {}
        for field in self.layer.fields():
            fieldsDict[field.name()] = field.typeName()
        return fieldsDict
    
    def bbox(self) -> list[float]:
        """
        Get the bounding box of the layer.
        
        It is of the form [xmin, ymin, xmax, ymax].
        """
        qgsRectangle = self.layer.extent()
        return [qgsRectangle.xMinimum(), qgsRectangle.yMinimum(), qgsRectangle.xMaximum(), qgsRectangle.yMaximum()]

    def size(self) -> int:
        """
        Get the number of features in the layer.
        """
        return self.layer.featureCount()
    
    def features(self, filter:str = None, bbox:list[float] = None, geometryfilter:HGeometry = None) -> list[HFeature]:
        """
        Get the features of the layer.

        The filter can be a SQL expression.
        The bbox is of the form [xmin, ymin, xmax, ymax].
        The geometry filter is a HGeometry object.

        Filter elements are exclusive, i.e. if filter is set, bbox and geometryfilter are ignored.
        """
        features = None
        if filter:
            features = self.layer.getFeatures(QgsFeatureRequest().setFilterExpression(filter))
        elif bbox:
            rect = QgsRectangle(*bbox)
            features = self.layer.getFeatures(QgsFeatureRequest().setFilterRect(rect))
        elif geometryfilter:
            features = self.layer.getFeatures(QgsFeatureRequest().setFilterRect(geometryfilter.geometry.boundingBox()))
            # then filter out manually by looping and intersecting
            features = [f for f in features if f.geometry().intersects(geometryfilter.get_qgs_geometry())]
        else:
            features = self.layer.getFeatures()

        return [HFeature(f) for f in features]

    def add_feature(self, geometry: HGeometry, attributes: list[any]):
        """
        Add a feature to the layer.
        """
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry(geometry.geometry.clone()))
        feature.setAttributes(attributes)
        self.layer.dataProvider().addFeature(feature)
    
    def add_features(self, features: list[HFeature]):
        """
        Add a list of features to the layer.
        """
        self.layer.dataProvider().addFeatures([f.feature for f in features])

    def subset_filter(self, filter: str):
        """
        Set a subset filter for the layer.
        """
        self.layer.setSubsetString(filter)

    def dump_to_gpkg(self, path: str, overwrite: bool = False, encoding: str = "UTF-8"):
        """
        Dump the layer to a GeoPackage file.
        """
        saveOptions = QgsVectorFileWriter.SaveVectorOptions()
        saveOptions.driverName = 'GPKG'
        saveOptions.layerName = self.layer.name()
        if overwrite:
            saveOptions.actionOnExistingFile = QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile
        else:
            saveOptions.actionOnExistingFile = QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer
        saveOptions.fileEncoding = encoding
        transform_context = QgsProject.instance().transformContext()
        error = QgsVectorFileWriter.writeAsVectorFormatV3(self.layer, path, transform_context, saveOptions)
        if error[0] != QgsVectorFileWriter.NoError:
            return error[1]
        return None

    def dump_to_shp(self, path: str, encoding: str = "UTF-8"):
        """
        Dump the layer to a Shapefile.
        """
        saveOptions = QgsVectorFileWriter.SaveVectorOptions()
        saveOptions.driverName = 'ESRI Shapefile'
        saveOptions.fileEncoding = encoding
        transform_context = QgsProject.instance().transformContext()
        error = QgsVectorFileWriter.writeAsVectorFormatV3(self.layer, path, transform_context, saveOptions)
        if error[0] != QgsVectorFileWriter.NoError:
            return error[1]
        return None

    def set_style(self, style: HStyle):
        """
        Set the style of the layer.
        """
        if style.type == "line":
            properties = {
                'line_color': style['stroke_color'],
                'line_width': style['stroke_width'],
                'capstyle': 'round',
                'joinstyle': 'round'
            }
            symbol = QgsLineSymbol.createSimple(properties)
            self.layer.renderer().setSymbol(symbol)
        elif style.type == "point":
            properties = {
                'name': style['name'],
                'size': style['size'],
                'angle': style['angle'],
                'color': style['fill_color'],
                'outline_color': style['stroke_color'],
                'outline_width': style['stroke_width']
            }
            symbol = QgsMarkerSymbol.createSimple(properties)
            self.layer.renderer().setSymbol(symbol)
        elif style.type == "polygon":
            properties = {
                'color': style['fill_color'],
                'outline_color': style['stroke_color'],
                'outline_width': style['stroke_width'],
                'capstyle': 'round',
                'joinstyle': 'round'
            }
            symbol = QgsFillSymbol.createSimple(properties)
            self.layer.renderer().setSymbol(symbol)
        else:
            raise ValueError(f"Unsupported style type: {style.type}")
        self.layer.triggerRepaint()

    
class HMap:
    """
    Class to manage the map view.
    """

    @staticmethod
    def add_layer(layer: HVectorLayer):
        """
        Add a layer to the current active map.
        """
        QgsProject.instance().addMapLayer(layer.layer)

    @staticmethod
    def remove_layer(layer: HVectorLayer):
        """
        Remove a layer from the current active map.
        """
        QgsProject.instance().removeMapLayer(layer.layer.id())

    @staticmethod
    def remove_layer_by_name(name: str):
        """
        Remove a layer from the current active map by name.
        """
        layers = QgsProject.instance().mapLayersByName(name)
        for layer in layers:
            QgsProject.instance().removeMapLayer(layer.id())

    @staticmethod
    def layers_by_name(name: str) -> list[HVectorLayer]:
        """
        Get layers by name.
        """
        layers = []
        for layer in QgsProject.instance().mapLayersByName(name):
            layers.append(HVectorLayer(layer, True))
        return layers

        