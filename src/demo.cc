#include <iomanip> 
#include <iostream>
#include <cmath> 
#include <getopt.h>
#include "rpy.h"
#include "dashmm/dashmm.h"

// Here we create the evaluator object needed in this demo, which must be
// instantiated before the call to dashmm::init to allow registeration of the
// relevant actions with the runtime system. 
dashmm::Evaluator<dashmm::Bead, dashmm::Bead, dashmm::RPY, 
                  dashmm::FMM97> rpy_fmm97{}; 
dashmm::Evaluator<dashmm::Bead, dashmm::Bead, dashmm::RPY, 
                  dashmm::Direct> rpy_direct{}; 

// This type collects the input arguments to the program 
struct InputArguments {
  int num_beads{10000}; 
  double radius{0.001}; 
  std::string distribution{"cube"}; 
  int threshold{80}; 
  int accuracy{3}; 
  bool verify{false}; 
}; 

// Print usage information
void print_usage(char *program) {
  std::cout 
    << "Usage: " << program << " [options]\n"
    << "--num-beads=num                number of beads, default is 10000\n"
    << "--radius=num                   radius of bead, default is 0.001\n"
    << "--distribution=[cube|sphere]   distribution of beads, default is cube\n"
    << "--threshold=num                max beads per leaf box, default is 80\n"
    << "--accuracy=[3|6]               accurate digits, default is 3\n"
    << "--verify=[yes|no]              verify accuracy, default is false\n";
}

// Parse the command line arguments, overiding any defaults at the request of
// the user. 
bool read_arguments(int argc, char **argv, InputArguments &retval) {
  int opt = 0; 
  static struct option long_options[] = {
    {"num-beads", required_argument, 0, 'n'}, 
    {"radius", required_argument, 0, 'r'}, 
    {"distribution", required_argument, 0, 'd'}, 
    {"threshold", required_argument, 0, 's'}, 
    {"accuracy", required_argument, 0, 'a'}, 
    {"verify", required_argument, 0, 'v'}, 
    {"help", no_argument, 0, 'h'}, 
    {0, 0, 0, 0}
  }; 

  int long_index = 0; 
  while ((opt = getopt_long(argc, argv, "n:r:d:s:a:v:h",
                            long_options, &long_index)) != -1) {
    std::string verifyarg{};     
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
    case 'v':
      verifyarg = optarg; 
      retval.verify = (verifyarg == std::string{"yes"}); 
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

  if (retval.radius < 0.0) {
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
    std::cout 
      << "DASHMM accelerated Rotne-Prager-Yamakawa tensor interaction\n"
      << "\nParameters:\n"
      << std::setw(70) << std::left << "Number of beads:"
      << std::setw(10) << std::right << retval.num_beads << "\n"
      << std::setw(70) << std::left << "Radius of bead:"
      << std::setw(10) << std::right << retval.radius << "\n"
      << std::setw(70) << std::left << "Distribution of beads:"
      << std::setw(10) << std::right << retval.distribution << "\n"
      << std::setw(70) << std::left << "Refinement limit:" 
      << std::setw(10) << std::right << retval.threshold << "\n"
      << std::setw(70) << std::left << "Accurate digits:"
      << std::setw(10) << std::right << retval.accuracy << "\n";
  }

  return true; 
}

// Prepare the bead data and put it in the global address space 
dashmm::Array<dashmm::Bead> prepare_beads(const InputArguments &args) {
  dashmm::Bead *beads{nullptr}; 

  // Dole out the beads equally 
  int num_ranks = dashmm::get_num_ranks(); 
  int myrank = dashmm::get_my_rank(); 
  int num_beads = args.num_beads / num_ranks; 
  int remainder = args.num_beads % num_ranks; 
  int offset = 0; 
  
  if (myrank < remainder) {
    num_beads++; 
    offset = myrank * num_beads;
  } else {
    offset = myrank * num_beads + remainder; 
  }

  beads = new dashmm::Bead[num_beads]; 
 
  if (args.distribution == std::string{"cube"}) {
    // Cube 
    for (int i = 0; i < num_beads; ++i) {
      beads[i].index = offset + i; 
      double temp[6]; 
      for (int j = 0; j < 6; ++j) 
        temp[j] = (double) rand() / RAND_MAX - 0.5;
      beads[i].q[0] = temp[0]; 
      beads[i].q[1] = temp[1]; 
      beads[i].q[2] = temp[2]; 
      beads[i].position = dashmm::Point{temp[3], temp[4], temp[5]}; 
    }
  } else {
    // Sphere 
    for (int i = 0; i < num_beads; ++i) {
      beads[i].index = offset + i; 
      double temp[5]; 
      for (int j = 0; j < 5; ++j) 
        temp[j] = (double) rand() / RAND_MAX - 0.5;         
      beads[i].q[0] = temp[0];
      beads[i].q[1] = temp[1];
      beads[i].q[2] = temp[2];

      double stheta = sin(temp[3]) * M_PI; 
      double ctheta = cos(temp[3]) * M_PI; 
      double sphi = sin(temp[4]) * M_PI; 
      double cphi = cos(temp[4]) * M_PI; 
      beads[i].position = 
        dashmm::Point{stheta * sphi, stheta * cphi, ctheta}; 
    }
  }

  dashmm::Array<dashmm::Bead> retval{}; 
  auto err = retval.allocate(num_beads, beads); 
  assert(err == dashmm::kSuccess); 

  return retval; 
}

// Compute an error characteristic for the values computed with a multipole
// method, and the values computed with direct summation 
void compare_results(dashmm::Bead *beads, int num_beads, 
                     dashmm::Bead *exacts, int num_exacts) {
  if (dashmm::get_my_rank()) 
    return;

  // beads have been sorted by the index field, exacts have not 
  
  // Loop over the exact results and compare 
  double numerator = 0.0; 
  double denominator = 0.0; 
  for (int i = 0; i < num_exacts; ++i) {
    int j = exacts[i].index;
    numerator += pow(beads[j].value[0] - exacts[i].value[0], 2) + 
      pow(beads[j].value[1] - exacts[i].value[1], 2) + 
      pow(beads[j].value[2] - exacts[i].value[2], 2); 
    denominator += pow(exacts[i].value[0], 2) + 
      pow(exacts[i].value[1], 2) + pow(exacts[i].value[2], 2); 
  }

  double error = sqrt(numerator / denominator); 

  std::cout << "\nError\n" 
            << std::setw(70) << std::left << "L2 norm error:" 
            << std::right << std::scientific << std::setprecision(4) 
            << error << "\n"; 
}

// Main driver of the evaluation 
void perform_evaluation_test(const InputArguments &args) {
  dashmm::Array<dashmm::Bead> bead_handle = prepare_beads(args); 

  dashmm::FMM97<dashmm::Bead, dashmm::Bead, dashmm::RPY> fmm97{}; 
  
  // One could specify kernel parameters with four entries, including
  // (1) radius of the bead 
  // (2) Boltzmann constant 
  // (3) Absolute temperature
  // (4) Solvent viscosity 
  std::vector<double> kparams(1, args.radius); 

  std::cout << "\nFMM Evaluation\n"; 
  auto err = rpy_fmm97.evaluate(bead_handle, bead_handle, args.threshold, 
                                &fmm97, args.accuracy, &kparams); 
  assert(err == dashmm::kSuccess); 

  if (args.verify) {
    // Select a few beads for the direct comparison
    int test_count{0}; 

    if (dashmm::get_my_rank() == 0) {
      test_count = std::min(400, args.num_beads); 
    }
      
    // Get the results from the global address space
    auto beads = bead_handle.collect(); 
    if (beads) {
      std::sort(&beads[0], &beads[args.num_beads], 
                [](const dashmm::Bead &a, const dashmm::Bead &b) -> bool {
                  return (a.index < b.index); 
                });
    }

    // Copy the test beads into test_beads
    dashmm::Bead *test_beads{nullptr}; 
    if (test_count) {
      test_beads = new dashmm::Bead[test_count]; 

      for (int i = 0; i < test_count; ++i) {
        int idx = i * (args.num_beads / test_count); 
        test_beads[i] = beads[idx]; 
        test_beads[i].value[0] = 0.0; 
        test_beads[i].value[1] = 0.0; 
        test_beads[i].value[2] = 0.0;
      }
    }

    // Create array for test beads
    dashmm::Array<dashmm::Bead> test_handle{}; 
    err = test_handle.allocate(test_count, test_beads); 
    assert(err == dashmm::kSuccess); 

    // Do direct evaluation 
    std::cout << "\nDirect Evaluation:\n";
    dashmm::Direct<dashmm::Bead, dashmm::Bead, dashmm::RPY> direct{}; 
    err = rpy_direct.evaluate(bead_handle, test_handle, args.threshold, 
                              &direct, args.accuracy, &kparams);
    assert(err == dashmm::kSuccess);

    // Retrieve results 
    auto test_results = test_handle.collect(); 

    // Test error 
    compare_results(beads.get(), args.num_beads, 
                    test_results.get(), test_count); 

    err = test_handle.destroy(); 
    assert(err == dashmm::kSuccess); 
  }

  // Free up resources 
  assert(bead_handle.destroy() == dashmm::kSuccess); 
}
 
// Program entry point
int main(int argc, char **argv) {
  auto err = dashmm::init(&argc, &argv); 
  assert(err == dashmm::kSuccess); 

  // Parse the command line arguments
  InputArguments inputargs; 
  if (read_arguments(argc, argv, inputargs)) {
    perform_evaluation_test(inputargs); 
  }

  err = dashmm::finalize(); 
  assert(err == dashmm::kSuccess); 

  return 0; 
}
