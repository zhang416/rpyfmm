#include <iomanip>
#include <iostream>
#include <getopt.h>
#include "rpy.h"
#include "dashmm/dashmm.h"

// Here we create the evaluator object needed in this demo, which must be
// instantiated before the call to dashmm::init to allow registeration of the
// relevant actions with the runtime system. 
dashmm::Evaluator<dashmm::Bead, dashmm::Bead, dashmm::RPY, dashmm::FMM97> rpy{}; 

/*
// This type collects the input arguments to the program 
struct InputArguments {
  int num_beads{10000}; 
  double radius{0.001}; 
  std::string distribution{"cube"}; 
  int threshold{80}; 
  int accuracy{3}; 
}; 

// Print usage information
void print_usage(char *program) {
  std::cout << "Usage: " << program << "\n"
            << " --num-beads=num\n"
            << "   number of beads, default is 10000\n"
            << " --radius=num\n"
            << "   radius of each bead, default is 0.001\n"
            << " --distribution=[cube|sphere]\n"
            << "   distribution of the beads, default is cube\n"
            << " --threshold=num\n"
            << "   maximum number of beads in a leaf box, default is 80\n"
            << " --accuracy=[3|6]\n"
            << "   number of accurate digits of the results, default is 3\n";
}

// Parse the command line arguments, overiding any defaults at the request of
// the user. 
bool read_arguments(int argc, char **argv, InputArgument &retval) {
  int opt = 0; 
  static struct option long_options[] = {
    {"num-beads", required_argument, 0, 'n'}, 
    {"radius", required_argument, 0, 'r'}, 
    {"distribution", required_argument, 0, 'd'}, 
    {"threshold", required_argument, 0, 's'}, 
    {"accuracy", required_argument, 0, 'a'}, 
    {0, 0, 0, 0}
  }; 

  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, "n:r:d:s:a:h",
                            long_options, &long_index)) != -1) {
    switch (opt) {
    case 'n':
      retval.num_beads = atoi(optarg);
      break;
    case 'r':
      retval.radius = atof(optarg);
      break;
    case 'd':
      retval.distribution = optarg;
      break;
    case 's':
      retval.threshold = atoi(optarg);
      break;
    case 'a':
      retval.accuracy = atoi(optarg);
      break;
    case 'h':
      print_usage(argv[0]);
      return false;
    case '?':
      return false;
    }
  }

  // Test the inputs 
  if (retval.num_beads < 1) {
    fprintf(stderr, "Usage ERROR: number of beads must be positive.\n");
    return false;
  } 

  if (retva.radius < 0.0) {
    fprintf(stderr, "Usage ERROR: radius of bead must be positive.\n");
    return false;
  }

  if (retval.distribution != "cube" && retval.distribution != "sphere") {
    fprintf(stderr, "Usage ERROR: unknown distribution type.\n");
    return false;
  }

  if (retval.threshold < 0) {
    fprintf(stderr, "Usage ERROR: threshold must be positive.\n");
    return false;
  }

  if (retval.accuracy != 3 && retval.accuracy != 6) {
    fprintf(stderr, "Usage ERROR: unsupported accuracy input.\n");
    return false;
  }

  // Print out summary 
  if (dashmm::get_my_rank() == 0) {
    std::cout << "DASHMM accelerated Rotne-Prager-Yamakawa tensor interaction"
              << "\n\nParameters:\n"
              << std::setw(50) << "... Number of beads:" 
              << std::setw(30) << std::right << retval.num_beads << "\n"
              << std::setw(50) << "... Radius of bead:"
              << std::setw(30) << std::right << std::setprecision(5) 
              << std::scientific << retval.radius << "\n"
              << std::setw(50) << "... Distribution of beads:"
              << std::setw(30) << std::right << retval.distribution << "\n"
              << std::setw(50) << "... Refinement limit:"
              << std::setw(30) << std::right << retval.threshold << "\n"
              << std::setw(50) << "... Number of accurate digits:"
              << std::setw(30) << std::right << retval.accuracy << "\n";
  }
}
  
*/

int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv); 
  assert(err == dashmm::kSuccess); 

  /*
  // Parse the command line arguments
  InputArguments inputargs; 
  if (read_arguments(argc, argv, inputargs)) {
    // do 
  }
  */

  err = dashmm::finalize(); 
  assert(err == dashmm::kSuccess); 

  return 0; 
}
