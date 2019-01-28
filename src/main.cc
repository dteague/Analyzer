#include "Analyzer.h"
#include <csignal>
#include <sys/stat.h>
#define Q(x) #x
#define QUOTE(x) Q(x)
#include QUOTE(MYANA)
#include <sys/types.h>
#include <dirent.h>
 


bool do_break;
void KeyboardInterrupt_endJob(int signum) {
    do_break = true;
}

void usage() {
  std::cout << "./Analyzer infile.root outfile.root\n";
  std::cout << "or\n";
  std::cout << "./Analyzer -out outfile.root -in infile.root infile.root infile.root ...\n";
  std::cout << "or\n";
  std::cout << "./Analyzer -out outfile.root -in infile.root\n";
  std::cout << "Available options are:\n";
  std::cout << "-C: use a different config folder than the default 'PartDet'\n";
  std::cout << "-meta: add metadata to file\n";
  std::cout << "-t: run over 100 events\n";
  std::cout << "\n";

  exit(EXIT_FAILURE);
}

void parseCommandLine(int argc, char *argv[], std::vector<std::string> &inputnames, std::string &outputname, bool &testRun, std::string &configFolder, bool& useMeta) {
  if(argc < 3) {
    std::cout << std::endl;
    std::cout << "You have entered too little arguments, please type:\n";
    usage();
  }
  for (int arg=1; arg<argc; arg++) {
    //// extra arg++ are there to move past flags
    if (strcmp(argv[arg], "-t") == 0) {
      testRun = true;
      continue;
    }else if (strcmp(argv[arg], "-C") == 0) {
      configFolder=argv[arg+1];
      std::cout << "Analyser: ConfigFolder " << configFolder << std::endl;
      arg++;
      continue;
    }else if (strcmp(argv[arg], "-in") == 0) {
      arg++;
      while( arg<argc and (argv[arg][0] != '-')){
        inputnames.push_back(argv[arg]);
        std::cout << "Analyser: Inputfilelist " << inputnames.back() << std::endl;
        arg++;
      }
      arg--; /// to counteract arg++ that is in the for loop
      continue;
    }else if (strcmp(argv[arg], "-out") == 0) {
      outputname=argv[arg+1];
      std::cout << "Analyser: Outputfile " << outputname << std::endl;
      arg++;
      continue;
    } else if(strcmp(argv[arg], "-meta") == 0) {
      useMeta = true;
      arg++;
    } else if(argv[arg][0] == '-') {
      std::cout << std::endl;
      std::cout << "You entered an option that doesn't exist.  Please use one of the options:" << std::endl;
      usage();
    }else if(inputnames.size()==0){
      inputnames.push_back(argv[arg]);
      std::cout << "Analyser: Inputfilelist " << inputnames.back() << std::endl;
    }else if(outputname==""){
      outputname = argv[arg];
      std::cout << "Analyser: Outputfile " << outputname << std::endl;
    }
  }

  if(inputnames.size() == 0) {
    std::cout << std::endl;
    std::cout << "No input files given!  Please type:" << std::endl;
    usage();
  } else if(outputname == "") {
    std::cout << std::endl;
    std::cout << "No output file given!  Please type:" << std::endl;
    usage();
  }


  
  //for( auto file: inputnames) {
    //ifstream ifile(file);
    //if ( !ifile && file.find("root://") == std::string::npos && file.find("root\\://") == std::string::npos) {
      //std::cout << "The file '" << inputnames.back() << "' doesn't exist" << std::endl;
      //exit(EXIT_FAILURE);
    //}
  //}
  return;
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

std::vector<std::string> read_directory(const std::string& name) {
  std::vector<std::string> v;
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != NULL) {
    v.push_back(dp->d_name);
  }
  closedir(dirp);
  return v;
}

std::vector<std::string> inputDirectory(std::vector<std::string> oldInput) {
  std::vector<std::string> newInput;
  for(auto filename: oldInput) {
    if (is_dir(filename.c_str())) {
      for(auto insideName : read_directory(filename)) {
	if(insideName.find(".root") != std::string::npos) 
	  newInput.push_back(filename+ "/" + insideName);
      }
    } else if(filename.find(".root") != std::string::npos) 
      newInput.push_back(filename);
  }
  return newInput;
}



int main (int argc, char* argv[]) {

  bool testRun = false;
  do_break =false;

  std::string outputname;
  std::string configFolder="PartDet";
  std::vector<std::string> inputnames;
  bool useMeta = false;

  //get the command line options in a nice loop
  parseCommandLine(argc, argv, inputnames, outputname, testRun, configFolder, useMeta);
  inputnames = inputDirectory(inputnames);

  //setup the analyser
  Analyzer testing(inputnames, outputname, configFolder);
  // SpechialAnalysis spechialAna = SpechialAnalysis(&testing);
  // spechialAna.init();

  //catch ctrl+c and just exit the loop
  //this way we still have the output
  signal(SIGINT,KeyboardInterrupt_endJob);

  size_t Nentries=testing.nentries;
  if(testRun){
    Nentries=100;
    testing.nentries=100;
  }
  //main event loop
  for(size_t i=0; i < Nentries; i++) {
    // if(i==0){
    //   spechialAna.begin_run();
    // }
    testing.clear_values();
    //    std::cout << "    testing.preprocess(i);" << std::endl;
    bool valid = testing.preprocess(i);
    if( !valid) continue;
    //    std::cout << "    testing.fill_histogram();" << std::endl;
    testing.fill_histogram();
    //    std::cout << "    spechialAna.analyze();" << std::endl;
    // spechialAna.analyze();
    //this will be set if ctrl+c is pressed
    if(do_break){
      testing.nentries=i+1;
      break;
    }
  }
  testing.printCuts();
  // spechialAna.end_run();
  return 0;
}
