
// METU Department of Mechanical Engineering
// ME 705 CFD for Incompressible Flows
// 2D INS solver based on unsteady vorticity - stream function formulation.
// It is solving the lid driven cavity (LDC) problem.
// Uses forward time and central space (FTCS) scheme.

// Dr. Cuneyt Sert


#include <iostream>
#include <cmath>
#include <fstream>    // To read/write from/to a file
#include <chrono>     // To calculate elapsed time
#include <iomanip>    // To control output precision


int main()
{
  auto time_begin = std::chrono::high_resolution_clock::now();

  double L = 1.0;
  double Ulid = 1.0;

  int nx = 51;
  int ny = nx;

  double dt = 0.005;
  int nStep = 50000;

  double visc = 1.0/1000.0;

  int maxGSiter = 3;
  double GStol = 1e-3;
  double steadyStateTol = 1e-5;
 
  double psi [nx][ny];
  double omega [nx][ny];
  double psiOld [nx][ny];
  double omegaOld [nx][ny];
  double u [nx][ny];
  double v [nx][ny];
  double residual[nx][ny];

  double h = L/(nx-1);
  double t = 0.0;

  double maxOmegaDiff;
  double maxPsiDiff;


  // We want nx to be even
  if (nx % 2 == 0) {
    std::cout << "\nError: nx needs to be odd! \n\n";
    return 1;
  }

  // Initialization
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
       psi[i][j] = 0;
       omega[i][j] = 0;
       psiOld[i][j] = 0;
       omegaOld[i][j] = 0;
       u[i][j] = 0;
       v[i][j] = 0;
    }
  }


  // Set the u at the top boundary once and don't change it later.
  for (int i = 0; i < nx; i++) {
    u[i][ny-1] = Ulid;
  }



  // Time loop
  for (int n = 1; n <= nStep; n++) {

    int GSiter;

    // Solve the Poison equation for stream function using Gauss-Seidel
    for (GSiter = 1; GSiter <= maxGSiter; GSiter++) {
      for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
          psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + h*h*omega[i][j]);
        }
      }

      // Calculate the GS residual
      for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
          residual[i][j] = psi[i][j] - 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] + h*h * omega[i][j]);
        }
      }

      // Use the L2 norm of the residuals for convergence check
      double sum = 0;
      for (int i = 1; i < nx-1; i++) {
        for (int j = 1; j < ny-1; j++) {
          sum = sum + residual[i][j]*residual[i][j];
        }
      }
      sum = std::sqrt(sum);

      if (sum < GStol) {
        break;
      }

    }  // End of Gauss-Seidel loop


    // Calculate vorticity at the boundaries
    for (int i = 0; i < nx; i++) {
      omega[i][0]  = -2 * psi[i][1] / (h*h);                    // Bottom wall
      omega[i][ny-1] = -2 * psi[i][ny-2] / (h*h) - 2*Ulid/h;    // Top wall
    }

    for (int j = 0; j < ny; j++) {
      omega[0][j] = -2 * psi[1][j] / (h*h);                     // Left wall
      omega[nx-1][j] = -2 * psi[nx-2][j] / (h*h);               // Right wall
    }

    // Calculate velocity components at the internal nodes by differentiating psiOld.
    // These will be used in the vorticity transport equation.
    for (int i = 1; i < nx-1; i++) {
      for (int j = 1; j < ny-1; j++) {
        u[i][j] =  (psiOld[i][j+1] - psiOld[i][j-1]) / (2*h);
        v[i][j] = -(psiOld[i+1][j] - psiOld[i-1][j]) / (2*h);
      }
    }

    // Calculate new vorticity values at the inner nodes using the vorticity transport equation.
    for (int i = 1; i < nx-1; i++) {
      for (int j = 1; j < ny-1; j++) {
        omega[i][j] = omegaOld[i][j] + dt * (
                     -0.5/h * (u[i][j] * (omegaOld[i+1][j] - omegaOld[i-1][j]) +
                               v[i][j] * (omegaOld[i][j+1] - omegaOld[i][j-1]))
                     + visc/(h*h) * (omegaOld[i+1][j] + omegaOld[i-1][j] + omegaOld[i][j+1] +
                                     omegaOld[i][j-1] - 4*omegaOld[i][j]));
      }
    }


    // Print information on the screen
    std::cout << "n = " << n
              << ",    t = "            << std::fixed << std::setprecision(5) << t + dt
              << ",    GSiter = "       << std::min(GSiter, maxGSiter)
              << ",    omega center = " << std::fixed << std::setprecision(5) << omega[(nx-1)/2][(ny-1)/2]
              << "\n";


    // Compare old and new omega and psi values to see whether state state is
    // reached or not.

    maxOmegaDiff = 0;
    maxPsiDiff = 0;
    double omegaDiff;
    double psiDiff;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        omegaDiff = std::abs(omega[i][j] - omegaOld[i][j]);
        if (omegaDiff > maxOmegaDiff) {
          maxOmegaDiff = omegaDiff;
        }

        psiDiff = std::abs(psi[i][j] - psiOld[i][j]);
        if (psiDiff > maxPsiDiff) {
          maxPsiDiff = psiDiff;
        }
      }
    }

    if (maxOmegaDiff <= steadyStateTol && maxPsiDiff <= steadyStateTol) {
      std::cout << "\nSteady state is reached. Exiting the time loop...\n";
      break;
    }


    // Get ready for the new time level
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        psiOld[i][j] = psi[i][j];
        omegaOld[i][j] = omega[i][j];
      }
    }

    t = t + dt;

  }  // End of the time loop


  // Write out the elapsed time
  auto time_end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_begin).count();
  std::cout << "Elapsed Time: " << duration << " [ms] \n\n";


  // Write out the variation of u at the vertical centerline of the cavity
  std::fstream output_file;
  output_file.open("LDC_output.txt", std::ios::out);
  output_file << "# LDC solution \n";
  output_file << "# nx      = " << nx   << "\n";
  output_file << "# dt      = " << dt   << "\n";
  output_file << "# visc    = " << visc << "\n";
  output_file << "# n       = " << t/dt  << "\n";
  output_file << "# t_final = " << t    << "\n";
  output_file << "#  y           u \n";

  for (int j = 0; j < ny; j++) {
    output_file <<           std::fixed << std::setprecision(5) << j*h 
                << "    " << std::fixed << std::setprecision(5) << u[(nx-1)/2][j] << "\n";
  }
  output_file.close();



  std::cout << "\nEnd of the main function is reached. Stopping.\n\n";
  return 0;
}
