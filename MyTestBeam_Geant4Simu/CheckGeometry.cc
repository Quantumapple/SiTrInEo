void CheckGeometry(const char* file_gdml = "config/run26999_LadderParallelToBieldAxis.cfg.gdml") {

  TCanvas* cgeometry3D = new TCanvas("cgeometryMC3D", "Telescope geometry 3D (gdml)");
  cgeometry3D->cd();
  TGeoManager::Import(file_gdml);
  gGeoManager->GetTopVolume()->Draw("ogl");
  
  return;
  
}