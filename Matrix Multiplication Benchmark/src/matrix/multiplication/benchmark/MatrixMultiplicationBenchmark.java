/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matrix.multiplication.benchmark;

/**
 *
 * @author Terry Wong
 */
public class MatrixMultiplicationBenchmark {

    private static final int NMIN = 50;
    private static final int NMAX = 500;
    private static final int STEP = 50;
    private static final int PERMUTATIONS = 48;
    private static double[][] a, b, c;

// Find the mininum of runtimes
    double min(double RunTime[]) {
        double min = RunTime[0];
        int n = PERMUTATIONS;

        for (int i = 1; i < n; i++) {
            if (RunTime[i] < min) {
                min = RunTime[i];
            }
        }

        return min;
    }

// Find the maximum of runtimes
    double max(double RunTime[]) {
        double max = RunTime[0];
        int n = PERMUTATIONS;

        for (int i = 1; i < n; i++) {
            if (RunTime[i] > max) {
                max = RunTime[i];
            }
        }

        return max;
    }

// Find the geometric mean of R(n)s
    double geometric_mean(double R[], int n) {
        double Rave, product = R[0];

        for (int i = 1; i < n; i++) {
            product *= R[i];
        }
        Rave = Math.pow((product), (1. / n));

        return Rave;
    }

    void allocateMatrix(int n) {
        // Allocate memory for matrices
        a = new double[n][n];
        b = new double[n][n];
        c = new double[n][n];
    }

    void deleteMatrix(int n) {
        // De-Allocate memory to prevent memory leak
        a = null;
        b = null;
        c = null;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        int i, j, k, M, m, n, nmin, nmax, step, row, index;
        double starttime, stoptime, Tmin, Tmax, R[], Q, Rave, Qave, v;
        double Tijk[] = new double[8], Tjik[] = new double[8], Tikj[] = new double[8],
                Tkij[] = new double[8], Tjki[] = new double[8], Tkji[] = new double[8];
        double RUN_TIME[] = new double[PERMUTATIONS];

        nmin = NMIN;
        nmax = NMAX;
        step = STEP;
        System.out.print(String.format("nmin: %d\tnmax: %d\tstep: %d\n", nmin, nmax, step));

        MatrixMultiplicationBenchmark mmb = new MatrixMultiplicationBenchmark();
        n = nmin;
        // Initialize a[][] and b[][] // Typical size from 100*100 to 500*500
        mmb.allocateMatrix(n);

        // Define arbitrary initial values of matrices a[][] and b[][]
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a[i][j] = Math.random();
                b[i][j] = Math.random();
            }
        }

        R = new double[nmax / nmin];
        index = 0;
        for (n = nmin; n <= nmax; n += step) {
            System.out.print(String.format("==========================================================================================================\n"));
            System.out.print(String.format("n=%-7d\t%-7s\t\t%-7s\t\t%-7s\t\t%-7s\t\t%-7s\t\t%-7s\n", n, "Tijk", "Tjik", "Tikj", "Tkij", "Tjki", "Tji"));
            System.out.print(String.format("==========================================================================================================\n"));
            mmb.allocateMatrix(n);
            M = (nmax * nmax * nmax) / (n * n * n);
            for (i = 0; i < n; i++) // Matrix [][] initialization
            {
                for (j = 0; j < n; j++) {
                    c[i][j] = 0.0;
                }
            }
            row = 0;

            /* Permutation 1: c[i][j] += a[i][k]*b[k][j]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[i][k] * b[k][j]; // Remove comment to make
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 2: c[i][j] += a[i][k]*b[j][k]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 3: c[i][j] += a[k][i]*b[k][j]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 4: c[i][j] += a[k][i]*b[j][k]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[i][j] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 5: c[j][i] += a[i][k]*b[k][j]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 6: c[j][i] += a[i][k]*b[j][k]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[i][k] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 7: c[j][i] += a[k][i]*b[k][j]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[k][i] * b[k][j];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification
            row++;

            /* Permutation 8: c[j][i] += a[k][i]*b[j][k]; */
            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // Basic i-j-k form (first of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tijk[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-i-k form
                {
                    for (i = 0; i < n; i++) {
                        for (k = 0; k < n; k++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjik[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (i = 0; i < n; i++) // i-k-j form
                {
                    for (k = 0; k < n; k++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tikj[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // k-i-j form
                {
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkij[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (j = 0; j < n; j++) // j-k-i form
                {
                    for (k = 0; k < n; k++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tjki[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            starttime = System.currentTimeMillis();
            for (m = 0; m < M; m++) // Repeat M times
            {
                for (k = 0; k < n; k++) // Final k-j-i form (last of 6 permutations)
                {
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) {
                            c[j][i] += a[k][i] * b[j][k];
                        }
                    }
                }
            }
            stoptime = System.currentTimeMillis();
            Tkji[row] = (stoptime - starttime) / M;
            v = c[n - 1][n - 1]; // Verification

            for (i = 0; i < 8; i++) {
                System.out.print(String.format("perm: %-5d\t%-7f\t%-7f\t%-7f\t%-7f\t%-7f\t%-7f\n", i + 1, Tijk[i], Tjik[i], Tikj[i], Tkij[i], Tjki[i], Tkji[i]));
            }

            i = 0;
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tijk[j];
            }
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tjik[j];
            }
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tikj[j];
            }
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tkij[j];
            }
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tjki[j];
            }
            for (j = 0; j < 8; j++) {
                RUN_TIME[i++] = Tkji[j];
            }

            Tmin = mmb.min(RUN_TIME);
            Tmax = mmb.max(RUN_TIME);
            R[index] = Tmax / Tmin;
            Q = 100 * (R[index] - 1) / (R[index] + 1);
            System.out.print(String.format("\nTmin: %f ms\tTmax: %f ms\tR(n): %f\tQ(n): %f%%\n", Tmin, Tmax, R[index], Q));
            index++;
            mmb.deleteMatrix(n);
        }

        System.out.print(String.format("\n********************************************************************************************************\n"));
        Rave = mmb.geometric_mean(R, index);
        Qave = 100 * (Rave - 1) / (Rave + 1);
        System.out.print(String.format("\nRave: %f\tQave: %f%%\n", Rave, Qave));
        System.out.print(String.format("End of the Matrix Multiplication Benchmark program!\n"));
    }

}
