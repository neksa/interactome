Clazz.declarePackage ("J.shapespecial");
Clazz.load (null, "J.shapespecial.Polyhedron", ["JU.BS", "$.P3i", "JU.Escape"], function () {
c$ = Clazz.decorateAsClass (function () {
this.modelIndex = 0;
this.centralAtom = null;
this.vertices = null;
this.ptCenter = 0;
this.visible = false;
this.normixes = null;
this.planes = null;
this.visibilityFlags = 0;
this.collapsed = false;
this.faceCenterOffset = 0;
this.distanceFactor = 0;
this.isFullyLit = false;
this.isValid = true;
this.colixEdge = 0;
Clazz.instantialize (this, arguments);
}, J.shapespecial, "Polyhedron");
Clazz.makeConstructor (c$, 
function (centralAtom, ptCenter, nPoints, planeCount, otherAtoms, normixes, planes, collapsed, faceCenterOffset, distanceFactor) {
this.centralAtom = centralAtom;
this.modelIndex = centralAtom.mi;
this.ptCenter = ptCenter;
this.vertices =  new Array (nPoints);
this.visible = true;
this.normixes =  Clazz.newShortArray (planeCount, 0);
this.planes =  new Array (planeCount);
for (var i = nPoints; --i >= 0; ) this.vertices[i] = otherAtoms[i];

for (var i = planeCount; --i >= 0; ) this.normixes[i] = normixes[i];

for (var i = planeCount; --i >= 0; ) {
var p = planes[i];
this.planes[i] = JU.P3i.new3 (p.x, p.y, p.z);
}
this.collapsed = collapsed;
this.faceCenterOffset = faceCenterOffset;
this.distanceFactor = distanceFactor;
}, "JM.Atom,~N,~N,~N,~A,~A,~A,~B,~N,~N");
Clazz.defineMethod (c$, "getState", 
function () {
var bs =  new JU.BS ();
for (var i = 0; i < this.ptCenter; i++) bs.set ((this.vertices[i]).i);

return "  polyhedra ({" + this.centralAtom.i + "}) to " + JU.Escape.eBS (bs) + (this.collapsed ? " collapsed" : "") + " distanceFactor " + this.distanceFactor + " faceCenterOffset " + this.faceCenterOffset + (this.isFullyLit ? " fullyLit" : "") + ";" + (this.visible ? "" : "polyhedra off;") + "\n";
});
});
