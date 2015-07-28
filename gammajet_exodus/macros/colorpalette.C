void colorpalette() { // Example of creating new colors (purples)
const Int_t colNum = 10; // and defining of a new palette
Int_t palette[colNum];
for (Int_t i=0; i<colNum; i++) {
// get the color and if it does not exist create it
if (! gROOT->GetColor(230+i) ){
TColor *color = new TColor(230+i,1-(i/((colNum)*1.0)),0.3,0.5,"");
} else {
TColor *color = gROOT->GetColor(230+i);
color->SetRGB(1-(i/((colNum)*1.0)),0.3,0.5);
}
palette[i] = 230+i;
}
gStyle->SetPalette(colNum,palette);
TF2 *f2 = new TF2("f2","exp(-(x^2)-(y^2))",-3,3,-3,3);
// two contours less than the number of colors in palette
f2->SetContour(colNum-2);
f2->Draw("cont");
}
