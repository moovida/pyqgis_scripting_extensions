# pyQGIS Scripting Extension

## Introduction

pyQGIS in a great and powerful tool to automate GIS tasks. 

But more than a scripting language it is a replication in python of the QGIS API. 
This is amazing but can be overwhelming for beginners.

This extension provides a set of tools to make scripting easier and userfriendly.

For example the following pyQGIS code to style a layer and set labels:

```python
# marker style
properties = {
    'name': 'square',
    'size': 8,
    'color': '255,0,0,128',
    'outline_color': 'black',
    'outline_width': 1,
    'angle': 45
}
symbol = QgsMarkerSymbol.createSimple(properties) 4
layer.renderer().setSymbol(symbol)

# label style
settings = QgsPalLayerSettings()
format = QgsTextFormat()
format.setFont(QFont('Arial', 8)) # font with size
format.setColor(QColor('red')) # font color
buffer = QgsTextBufferSettings()
buffer.setEnabled(True)
buffer.setSize(0.50)
buffer.setColor(QColor('black'))
# create a halo around the text
format.setBuffer(buffer) # set the halo in the text format
settings.setFormat(format) # set the text format in the layer settings
settings.fieldName = "NAME" # use the NAME field as label
# label positioning with offsets
settings.placement = QgsPalLayerSettings.OverPoint
settings.xOffset = 0.0
settings.yOffset = -8.5
# enable labels on the layer
labels = QgsVectorLayerSimpleLabeling(settings)
layer.setLabelsEnabled(True)
layer.setLabeling(labels)
```

boils down to this:

```python
from pyqgis_scripting_ext.core import *

pointStyle = HMarker("square", 8, 45) + HFill("0,255,0,128") + HStroke("black", 1) + HLabel("NAME", yoffset = -5) + HHalo("white", 1)
layer.set_style(pointStyle)
```

Also the API is a bit difficult to follow when it comes to the main geometry objects. There are QgsGeometry, QgsPointXY and the 
specializations of QgsAbstractGeometry like QgsPolygon, QgsLineString or QgsPoint. Here we use a single Geometry object with a 
common parent. 


And then there is this other thing. Since pyQGIS is bound to its lowlevel API, memory management can be an issue. 
If you reuse the same variable name for different objects, most of the times you will experience a QGIS crash.

This is quite brutal and a huge problem for beginners. Take the following example:

```python
p = QgsPoint(12, 46)
pg = QgsGeometry(p)

coords = [[31,11], [10,30], [20,40], [40,40]]
l = QgsLineString([ QgsPoint(p[0], p[1]) for p in coords])
lg = QgsGeometry(l)

print(pg)
print(lg)

# assume somewhere else in the code you just need to recreate 
# the QgsGeometry for whatever reason, using the previous p object
pg = QgsGeometry(p)
print(pg)
```

This super simple script crashes QGIS.

In this extensions an effort is made to avoid this issue by cloning/recreating objects when needed.


Clearly this library is limited against everything you can do in pyQGIS, but it is a starting point. 

One last thing: we chose simplicity over efficiency, so you might experience performance issues for large datasets.

## Compatibility

Compatibility table:

| QGIS Version | pyQGIS Scripting Extension |
|---------------|---------------------------|
| Prizren          | 0.1.0                     |

## Installation

At the time being the lib is available in pypi. You can install it with pip:

```bash
pip install pyqgis-scripting-ext
```

To make sure you install it using the python used by QGIS, you can run directly from the QGIS python console:

```python
import subprocess
subprocess.run(["pip", "install", "pyqgis-scripting-ext"])
```

