{
// Print ROOT version number
  Int_t ver, mver,plevel;
  sscanf(gROOT->GetVersion(),"%d.%d/%d",&ver,&mver,&plevel);
  printf("#define __ROOT_VERSION__ %d\n",ver);
  printf("#define __ROOT_MINORVERSION__ %d\n",mver);
  printf("#define __ROOT_PATCHLEVEL__ %d\n",plevel);
  Int_t fulvrs=10000*ver+100*mver+plevel;
  printf("#define __ROOT_FULLVERSION__ %d\n",fulvrs);
}
