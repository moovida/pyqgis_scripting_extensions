from qgis.core import *
from qgis.gui import *
import os
from PyQt5.QtGui import QColor, QFont

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
    def __init__(self, field:str = None, font: str = "Arial", color: str = 'black', size: float = 10.0,  \
                 xoffset: float = 0.0, yoffset: float = 0.0, along_line: bool = False, bold: bool = False, italic: bool = False):    
        self.font = font
        self.color = color
        self.size = size
        self.field = field
        self.xoffset = xoffset
        self.yoffset = yoffset
        self.along_line = along_line
        self.bold = bold
        self.italic = italic
        properties = {
            "label_font": font,
            "label_color": color,
            "label_size": size,
            "label_field": field,
            "label_xoffset": xoffset,
            "label_yoffset": yoffset,
            "label_along_line": along_line,
            "label_bold": bold,
            "label_italic": italic
        }
        super().__init__(properties)

class HHalo(HStyle):
    def __init__(self, color: str = 'white', width: float = 1):
        self.width = width
        self.color = color
        properties = {
            "halo_width": width,
            "halo_color": color,
        }
        super().__init__(properties)


class HFeature:
    def __init__(self, feature: QgsFeature):
        self.feature = feature
        self.geometry = HGeometry.from_specialized(self.feature.geometry())
        self.attributes = self.feature.attributes()

    # @staticmethod
    # def create(geometry: HGeometry, attributes: list[any]):
    #     """
    #     Create a new feature based on the given geometry and attributes.

    #     The attributes must be a dictionary with the field names as keys 
    #     and the same order as the schema they will be written to..
    #     """
    #     feature = QgsFeature()
    #     feature.setGeometry(QgsGeometry(geometry.geometry.clone()))
    #     feature.setAttributes(attributes)
    #     return HFeature(feature)
    
    def __str__(self):
        wkt = self.geometry.asWkt()
        maxLength = 50
        if len(wkt) > maxLength:
            wkt = wkt[:maxLength] + "..."
        string = f"Geometry: {wkt} \nAttributes: {self.attributes}"
        return string

class HLayer:
    def __init__(self, layer) -> None:
        self.layer = layer

class HRasterLayer(HLayer):
    def __init__(self, layer: QgsRasterLayer):
        super().__init__(layer)

    

class HVectorLayer(HLayer):
    """
    A wrapper for a vector layer.

    Allows to open or create new memory layers.
    """
    def __init__(self, layer: QgsVectorLayer, isReadOnly: bool):
        super().__init__(layer)
        self.isReadOnly = isReadOnly
        self.prjcode = self.layer.crs().authid()
        self.fields = {}
        self.field_names = []
        for field in self.layer.fields():
            self.fields[field.name()] = field.typeName()
            self.field_names.append(field.name())
        

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
    def new(name: str, geometry_type: str, prjcode: str, fields: dict[str, str]):
        """
        Create a new vector layer in memory.
        """
        definition = f"{geometry_type}?crs={prjcode}"
        for fname, ftype in fields.items():
            ftypel = ftype.lower()
            if ftypel.startswith("int"):
                ftype = "integer"
            elif ftypel.startswith("string"):
                ftype = "string"
            elif ftypel.startswith("float"):
                ftype = "float"
            elif ftypel.startswith("double") or ftypel.startswith("real"):
                ftype = "double"
            definition += f"&field={fname}:{ftype}"
        definition += "&index=yes"
        layer = QgsVectorLayer(definition, name, "memory")
        return HVectorLayer(layer, False)
    
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
    
    def field_index(self, field_name: str) -> int:
        """
        Get the index of a field by name.
        """
        return self.layer.fields().indexFromName(field_name)

    def add_feature(self, geometry: HGeometry, attributes: list[any]):
        """
        Add a feature to the layer.
        """
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry(geometry.geometry.clone()))
        feature.setAttributes(attributes)
        self.layer.dataProvider().addFeature(feature)

    def subset_filter(self, filter: str):
        """
        Set a subset filter for the layer.
        """
        self.layer.setSubsetString(filter)

    def sub_layer(self, intersection_geometry: HGeometry, name: str = None, fields: list[str] = None):
        """
        Create a sublayer from the intersection of the layer with the given geometry.
        """
        if name is None:
            name = self.layer.name()# + "_sublayer"
        features = self.features(geometryfilter=intersection_geometry)
        geomType = self.layer.geometryType()
        geomType = str(geomType).split('.')[-1]
        if geomType.startswith("Line"):
            geomType = "LineString"

        subLayerFields = self.fields
        if fields:
            subLayerFields = {f: self.fields[f] for f in fields}
        sublayer = HVectorLayer.new(name, geomType, self.prjcode, subLayerFields)

        fieldsIndexes = None
        if fields:
            fieldsIndexes = [self.field_index(f) for f in subLayerFields.keys()]

        # make a deep copy of the features
        for feature in features:
            if fieldsIndexes:
                attributes = [feature.attributes[i] for i in fieldsIndexes]
            else:
                attributes = feature.attributes.copy()
            sublayer.add_feature(feature.geometry, attributes)
            
        return sublayer

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
        renderer = self.layer.renderer()
        symbol = self.style_to_symbol(style)
        if renderer is None:
            renderer = QgsSingleSymbolRenderer(symbol)
            self.layer.setRenderer(renderer)
        renderer.setSymbol(symbol)

        self.set_labels(style)
        
        self.layer.triggerRepaint()

    def set_labels(self, style:HStyle):
        # if labels are set, add them
        if 'label_field' in style.properties:
            field = style.properties['label_field']
            font = "Arial"
            if 'label_font' in style.properties:
                font = style.properties['label_font']
            size = 14
            if 'label_size' in style.properties:
                size = style.properties['label_size']
            color = "black"
            if 'label_color' in style.properties:
                color = style.properties['label_color']
            bold = False
            if 'label_bold' in style.properties:
                bold = style.properties['label_bold']
            italic = False
            if 'label_italic' in style.properties:
                italic = style.properties['label_italic']

            settings = QgsPalLayerSettings()
            format = QgsTextFormat()
            format.setFont(QFont(font)) 
            format.setColor(QColor(color))
            format.setSize(size) # font size. The size in QFont is ignored
            format.setForcedBold(bold)
            format.setForcedItalic(italic)

            if 'halo_width' in style.properties:
                buffer = QgsTextBufferSettings()
                buffer.setEnabled(True)
                buffer.setSize(style.properties['halo_width'])
                buffer.setColor(QColor(style.properties['halo_color'] or "white"))
                format.setBuffer(buffer)

                # create a halo around the text
                format.setBuffer(buffer) # set the halo in the text format

            settings.setFormat(format) # set the text format in the layer settings
            settings.fieldName = field

            # make a simple check to see if brackets are contained. If they are, it is handled as an expression
            if "(" in field:
                settings.isExpression = True
            else:
                settings.isExpression = False

            # label positioning with offsets
            if style.type == "line" and 'label_along_line' in style.properties:
                settings.placement = QgsPalLayerSettings.Line
            elif style.type == "point": 
                settings.placement = QgsPalLayerSettings.OverPoint
                settings.xOffset = 0.0
                settings.yOffset = 0.0
                if 'label_xoffset' in style.properties:
                    settings.xOffset = style.properties['label_xoffset']
                if 'label_yoffset' in style.properties:
                    settings.yOffset = style.properties['label_yoffset']
            
            # enable labels on the layer
            labels = QgsVectorLayerSimpleLabeling(settings)
            self.layer.setLabelsEnabled(True)
            self.layer.setLabeling(labels)


    def style_to_symbol(self, style):
        symbol = None
        if style.type == "line":
            properties = {
                'line_color': style.properties.get('stroke_color') or "black",
                'line_width': style.properties.get('stroke_width') or 1,
                'capstyle': 'round',
                'joinstyle': 'round'
            }
            symbol = QgsLineSymbol.createSimple(properties)
        elif style.type == "point":
            properties = {
                'name': style.properties.get('marker_name') or "circle",
                'size': style.properties.get('marker_size') or 10,
                'angle': style.properties.get('marker_angle') or 0,
                'color': style.properties.get('fill_color') or "red"  ,
                'outline_color': style.properties.get('stroke_color') or "black",
                'outline_width': style.properties.get('stroke_width') or 1
            }
            symbol = QgsMarkerSymbol.createSimple(properties)
        elif style.type == "polygon":
            properties = {
                'color': style.properties.get('fill_color') or "red",
                'outline_color': style.properties.get('stroke_color') or "black",
                'outline_width': style.properties.get('stroke_width') or 1,
                'capstyle': 'round',
                'joinstyle': 'round'
            }
            symbol = QgsFillSymbol.createSimple(properties)
        else:
            raise ValueError(f"Unsupported style type: {style.type}")
        return symbol

    def set_graduated_style(self, field: str, ranges: list[tuple[float, float]], styles_list: list[HStyle], labelstyle: HStyle = None):
        """
        Set a graduated style to the layer.
        """
        renderer = QgsGraduatedSymbolRenderer(field)
    
        for (range, style) in zip(ranges, styles_list):
            min = range[0]
            max = range[1]
            label = ""
            if min != float('-inf'):
                label = f'{min} < '
            label += field
            if max != float('inf'):
                label += f' < {max}'
            if min == max:
                label = f'{field} = {min}'
            classificationRange = QgsClassificationRange(label, min, max)
            
            symbol = self.style_to_symbol(style)
            rendererRange = QgsRendererRange(classificationRange, symbol)
            renderer.addClassRange(rendererRange)

        if labelstyle:
            self.set_labels(labelstyle)

        self.layer.setRenderer(renderer)

        self.layer.triggerRepaint()

class HMapCanvas():
    def __init__(self, iface: QgisInterface = None):
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

    def set_layers(self, layers: list[HVectorLayer]):
        self.canvas.setLayers([layer.layer for layer in layers])

    def add_geometry(self, geometry: HGeometry, color: str = 'red', width: int = 2):
        if isinstance(geometry, HPoint):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Point)
            r.addGeometry(QgsGeometry(geometry.geometry.clone()))
            fillColor = QColor(color)
            fillColor.setAlphaF(0.5)
            r.setColor(fillColor)
            strokeColor = QColor(color)
            r.setStrokeColor(strokeColor)
            
            r.setWidth(width)
        elif isinstance(geometry, HMultiPoint):
            r = QgsRubberBand(self.canvas, Qgis.GeometryType.Point)
            for point in geometry.geometries():
                r.addGeometry(QgsGeometry(point.geometry.clone()))
            r.setFillColor(QColor(color))
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

    
class HMap:
    """
    Class to manage the map view.
    """

    @staticmethod
    def add_layer(layer: HLayer):
        """
        Add a layer to the current active map.
        """
        QgsProject.instance().addMapLayer(layer.layer)

    @staticmethod
    def remove_layer(layer: HLayer):
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
    def remove_layers_by_name(names: list[str]):
        """
        Remove a list of layers from the current active map by name.
        """
        for name in names:
            layers = QgsProject.instance().mapLayersByName(name)
            for layer in layers:
                QgsProject.instance().removeMapLayer(layer.id())

    @staticmethod
    def layers_by_name(name: str) -> list[HLayer]:
        """
        Get layers by name.
        """
        layers = []
        for layer in QgsProject.instance().mapLayersByName(name):
            layers.append(HLayer(layer))
        return layers
    
    @staticmethod
    def get_osm_layer() -> HRasterLayer:
        """
        Add an OpenStreetMap layer to the current active map.
        """
        return HRasterLayer(QgsRasterLayer("type=xyz&url=http://tile.openstreetmap.org/{z}/{x}/{y}.png", "OpenStreetMap", "wms"))
        
    
    @staticmethod
    def get_tms_layer(url: str, name: str) -> HRasterLayer:
        """
        Add a TMS layer to the current active map.
        """
        return HRasterLayer(QgsRasterLayer(f"type=xyz&url={url}", name, "wms"))



class HPrinter:
    
    def __init__(self, iface: QgisInterface, name:str = None, page:str = "A4") -> None:
        """
        Create a new printer object.

        If a name is given, it will be used to save the layout to the layout manager.
        """
        self.iface = iface
        self.project = QgsProject.instance()
        self.layout = QgsPrintLayout(self.project)
        self.layout.initializeDefaults() # add A4 page to layout by default

        # ls = QgsLayoutSize.A4
        # if page == "A3":
        #     ls = QgsLayoutSize.A3
        # elif page == "A5":
        #     ls = QgsLayoutSize.A5
        # layout_size = QgsLayoutSize(ls)
        # self.layout.pageCollection().pageLayout().setPageSize(layout_size)

        if name:
            # if you want to save the layout
            self.layout.setName(name)
            self.project.layoutManager().addLayout(self.layout)
        self.map = None
    
    def checkMap(self):
        if not self.map:
            raise ValueError("A map must be added to the layout first")

    def add_map(self, x:int = 0, y:int = 0, width:int = 296, height:int = 210, frame:bool = True, extent:list[float] = None):
        """
        Add a map to the layout.
        """
        self.map = QgsLayoutItemMap(self.layout)
        self.map.attemptMove(QgsLayoutPoint(x,y, QgsUnitTypes.LayoutMillimeters))
        self.map.attemptResize(QgsLayoutSize(width, height, QgsUnitTypes.LayoutMillimeters))
        self.map.setExtent(QgsRectangle(0, 0, width, height))

        ext = self.iface.mapCanvas().extent()
        if extent:
            ext = QgsRectangle(extent[0], extent[1], extent[2], extent[3])
        self.map.zoomToExtent(ext)
        self.map.setFrameEnabled(frame)
        self.layout.addLayoutItem(self.map)

    def add_legend(self, x:int = 230, y:int = 30, width:int = 150, height:int = 100, frame:bool = True, max_symbol_size:float = 3):
        """
        Add a legend to the layout.
        """
        self.checkMap()
        legend = QgsLayoutItemLegend(self.layout)
        legend.setLinkedMap(self.map) 
        legend.setFrameEnabled(frame)
        legend.attemptMove(QgsLayoutPoint(x,y, QgsUnitTypes.LayoutMillimeters))
        legend.attemptResize(QgsLayoutSize(width, height, QgsUnitTypes.LayoutMillimeters))
        legend.setMaximumSymbolSize(max_symbol_size)
        self.layout.addLayoutItem(legend)
    
    def add_scalebar(self, x:int = 10, y:int = 190, units:str = "km", style:str = "Single Box", segments:int = 2, \
                     unit_per_segment:int = 100, font:str = "Arial", font_size:int = 8):
        """
        Add a scalebar to the layout.
        """
        self.checkMap()
        scalebarTypes = ['Single Box', 'Double Box', 'Line Ticks Middle', 'Line Ticks Down', 'Line Ticks Up', 'Stepped Line', 'Hollow', 'Numeric']
        if not style in scalebarTypes:
            raise ValueError(f"Unsupported scalebar style: {style}. Supported types are: {scalebarTypes}")

        scalebar = QgsLayoutItemScaleBar(self.layout)
        scalebar.setStyle(style) #'Line Ticks Up')
        scalebar.setLinkedMap(self.map)

        f = QFont()
        f.setFamily(font)
        f.setPointSize(font_size)
        scalebar.setFont(f)
        # scalebar.applyDefaultSize()
        if units == "m":
            scalebar.setUnits(QgsUnitTypes.DistanceMeters)
        elif units == "km":
            scalebar.setUnits(QgsUnitTypes.DistanceKilometers)
        elif units == "mi":
            scalebar.setUnits(QgsUnitTypes.DistanceMiles)
        elif units == "ft":
            scalebar.setUnits(QgsUnitTypes.DistanceFeet)
        
        scalebar.setUnitLabel(units)
        scalebar.setNumberOfSegments(segments)
        scalebar.setUnitsPerSegment(unit_per_segment)
        scalebar.attemptMove(QgsLayoutPoint(x,y, QgsUnitTypes.LayoutMillimeters))
        self.layout.addLayoutItem(scalebar)

    def add_label(self, x:int = 10, y:int = 10, text:str = "a label", font:str = "Arial", font_size:int = 12, bold:bool = False, italic:bool = False, underline:bool = False):
        """
        Add a label to the layout.
        """
        label = QgsLayoutItemLabel(self.layout)
        label.setText(text)
        
        f = QFont()
        f.setFamily(font)
        f.setPointSize(font_size)
        if bold:
            f.setBold(True)
        if italic:
            f.setItalic(True)
        if underline:
            f.setUnderline(True)

        label.setFont(f)
        label.attemptMove(QgsLayoutPoint(x,y, QgsUnitTypes.LayoutMillimeters))
        self.layout.addLayoutItem(label)

    def dump_to_pdf(self, path:str):
        """
        Dump the layout to a PDF file.
        """
        exporter = QgsLayoutExporter(self.layout)
        settings = QgsLayoutExporter.PdfExportSettings()
        exporter.exportToPdf(path, settings)

    def dump_to_image(self, path:str, width:int = 800, height:int = 600, dpi:int = 300):
        """
        Dump the layout to an image file.
        """
        exporter = QgsLayoutExporter(self.layout)
        settings = QgsLayoutExporter.ImageExportSettings()
        settings.dpi = dpi
        settings.width = width
        settings.height = height
        exporter.exportToImage(path, settings)