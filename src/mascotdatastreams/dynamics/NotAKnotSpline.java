package mascotdatastreams.dynamics;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * Not-a-knot cubic spline interpolation implementation.
 * Adapted from CoalRe package for use in MascotDatastreams.
 * 
 * OVERVIEW:
 * This class provides cubic spline interpolation with not-a-knot boundary conditions,
 * which ensures smooth interpolation between data points. It also precomputes
 * values at grid points for efficient lookup.
 * 
 * CUBIC SPLINE MATHEMATICS:
 * A cubic spline consists of cubic polynomials S_i(t) on each interval [x[i], x[i+1]]:
 * S_i(t) = a[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])^2 + d[i]*(t-x[i])^3
 * 
 * COEFFICIENT INTERPRETATION:
 * - a[i] = y[i]: Function value at knot x[i] (given data)
 * - b[i]: Coefficient of linear term (t-x[i])
 * - c[i]: Coefficient of quadratic term (t-x[i])^2 
 * - d[i]: Coefficient of cubic term (t-x[i])^3
 * 
 * The actual derivatives at knots are computed from these coefficients:
 * - S'(x[i]) = b[i] (first derivative at left endpoint of interval i)
 * - S''(x[i]) = 2*c[i] (second derivative at left endpoint of interval i)
 * - S'''(x[i]) = 6*d[i] (third derivative at left endpoint of interval i)
 * 
 * NOT-A-KNOT CONDITION:
 * The not-a-knot condition is a boundary condition that eliminates the need
 * for additional constraints at the endpoints. It requires:
 * - At x[1]: S'''_0(x[1]) = S'''_1(x[1]) (third derivative continuity)
 * - At x[n-1]: S'''_{n-2}(x[n-1]) = S'''_{n-1}(x[n-1]) (third derivative continuity)
 * 
 * This creates a tridiagonal system of equations that is solved using
 * LU decomposition to determine the second derivatives at all knots.
 */
public class NotAKnotSpline {
    
    private final double[] x;  // knots (time points) - strictly increasing
    private final double[] y;  // function values at knots
    private double[][] splineCoeffs;  // spline coefficients [i][0]=d, [i][1]=c, [i][2]=b, [i][3]=a
    
    // NOTE: splineCoeffs[i] = [d, c, b, a] are coefficients of the cubic polynomial
    // For interval [x[i], x[i+1]], the polynomial is:
    // S_i(t) = a[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])^2 + d[i]*(t-x[i])^3
    // where splineCoeffs[i][3] = a[i], splineCoeffs[i][2] = b[i], etc.
    
    // Grid point functionality for efficient lookup
    private double[] gridPoints;
    private double[] gridValues;
    private double[] gridDerivatives;
    private double gridStep;
    private double gridStart;
    private double gridEnd;
    
    /**
     * Constructs a not-a-knot cubic spline from given data points.
     * 
     * @param x knots (must be strictly increasing)
     * @param y function values at knots
     * @throws IllegalArgumentException if x is not strictly increasing or arrays have different lengths
     */
    public NotAKnotSpline(double[] x, double[] y) {
        this(x, y, 1000); // Default to 1000 grid points
    }
    
    /**
     * Constructs a not-a-knot cubic spline from given data points with specified number of grid points.
     * 
     * @param x knots (must be strictly increasing)
     * @param y function values at knots
     * @param numGridPoints number of grid points for efficient lookup (must be >= 2)
     * @throws IllegalArgumentException if x is not strictly increasing, arrays have different lengths, or numGridPoints < 2
     */
    public NotAKnotSpline(double[] x, double[] y, int numGridPoints) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Arrays x and y must have the same length");
        }
        if (x.length < 2) {
            throw new IllegalArgumentException("At least 2 data points are required");
        }
        if (numGridPoints < 2) {
            throw new IllegalArgumentException("numGridPoints must be >= 2");
        }
        
        // Check that x is strictly increasing
        for (int i = 1; i < x.length; i++) {
            if (x[i] <= x[i-1]) {
                throw new IllegalArgumentException("x must be strictly increasing");
            }
        }
        
        this.x = x.clone();
        this.y = y.clone();
        
        int n = x.length;
        this.splineCoeffs = new double[n-1][4];
        
        computeCoefficients();
        setupGridPoints(numGridPoints);
    }
    
    /**
     * Computes the spline coefficients using proper not-a-knot boundary conditions.
     * This implementation uses the not-a-knot algorithm with matrix solving.
     * 
     * NOT-A-KNOT CONDITION:
     * The not-a-knot condition requires that the third derivative be continuous
     * at the second and second-to-last interior knots. This eliminates the need
     * for additional boundary conditions at the endpoints.
     * 
     * ALGORITHM:
     * 1. Set up tridiagonal system A*mu = r for not-a-knot conditions
     * 2. Solve for second derivatives mu using LU decomposition
     * 3. Compute polynomial coefficients from second derivatives
     * 
     * SPECIAL CASE FOR 3 POINTS:
     * When n=3, the not-a-knot condition creates a singular matrix.
     * In this case, we fall back to natural spline boundary conditions
     * (second derivative = 0 at endpoints).
     */
    private void computeCoefficients() {
        int n = x.length;
        
        // Handle edge case: only 2 points - linear interpolation
        if (n == 2) {
            double h = x[1] - x[0];
            double delta = (y[1] - y[0]) / h;
            splineCoeffs[0][0] = 0.0;  // d = 0 (no cubic term)
            splineCoeffs[0][1] = 0.0;  // c = 0 (no quadratic term)
            splineCoeffs[0][2] = delta; // b = slope
            splineCoeffs[0][3] = y[0];  // a = y[0]
            return;
        }

        // Calculate h values (difference between x values)
        double[] h = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = x[i + 1] - x[i];
        }

        // Calculate the slope between two knots
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        // Create the tridiagonal system
        RealMatrix A = new Array2DRowRealMatrix(n, n);
        RealVector r = new ArrayRealVector(n);

        if (n == 3) {
            // Special case: 3 points - use natural spline boundary conditions
            // Natural spline: second derivative = 0 at endpoints
            // This gives us mu[0] = 0 and mu[n-1] = 0
            
            // First row: natural boundary condition at x[0] (mu[0] = 0)
            A.setEntry(0, 0, 1.0);
            A.setEntry(0, 1, 0.0);
            A.setEntry(0, 2, 0.0);
            r.setEntry(0, 0.0);
            
            // Middle row: continuity condition at x[1]
            A.setEntry(1, 0, h[0]);
            A.setEntry(1, 1, 2 * (h[0] + h[1]));
            A.setEntry(1, 2, h[1]);
            r.setEntry(1, 6 * (delta[1] - delta[0]));
            
            // Last row: natural boundary condition at x[2] (mu[2] = 0)
            A.setEntry(2, 0, 0.0);
            A.setEntry(2, 1, 0.0);
            A.setEntry(2, 2, 1.0);
            r.setEntry(2, 0.0);
        } else {
            // General case: n >= 4 points - use not-a-knot boundary conditions
            // Set up the system A*mu = r for the not-a-knot condition
            // First row: not-a-knot condition at x[1]
            A.setEntry(0, 0, h[1]);
            A.setEntry(0, 1, -(h[0] + h[1]));
            A.setEntry(0, 2, h[0]);
            
            // Last row: not-a-knot condition at x[n-2]
            A.setEntry(n - 1, n - 3, h[n - 2]);
            A.setEntry(n - 1, n - 2, -(h[n - 2] + h[n - 3]));
            A.setEntry(n - 1, n - 1, h[n - 3]);

            for (int i = 1; i < n - 1; i++) {
                A.setEntry(i, i - 1, h[i - 1]);
                A.setEntry(i, i, 2 * (h[i - 1] + h[i]));
                A.setEntry(i, i + 1, h[i]);
                r.setEntry(i, 6 * (delta[i] - delta[i - 1]));
            }
        }

        // Solve for mu
        DecompositionSolver solver = new LUDecomposition(A).getSolver();
        RealVector mu = solver.solve(r);

        // Store coefficients for each spline segment
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i][0] = (mu.getEntry(i + 1) - mu.getEntry(i)) / (6 * h[i]);
            splineCoeffs[i][1] = mu.getEntry(i) / 2;
            splineCoeffs[i][2] = delta[i] - h[i] * (2 * mu.getEntry(i) + mu.getEntry(i + 1)) / 6;
            splineCoeffs[i][3] = y[i];
        }
    }
    
    /**
     * Sets up grid points for efficient value lookup.
     * Precomputes value at regular intervals.
     * 
     * @param numGridPoints number of grid points to create
     */
    private void setupGridPoints(int numGridPoints) {
        // Use exact knot range (matches CoalRe implementation)
        gridStart = x[0];
        gridEnd = x[x.length-1];
        
        gridStep = (gridEnd - gridStart) / (numGridPoints - 1);
        
        gridPoints = new double[numGridPoints];
        gridValues = new double[numGridPoints];
        gridDerivatives = new double[numGridPoints];
        
        // Precompute values and derivatives at grid points
        for (int i = 0; i < numGridPoints; i++) {
            double t = gridStart + i * gridStep;
            gridPoints[i] = t;
            gridValues[i] = value(t);
            gridDerivatives[i] = derivative(t); 
        }
    }
    
    /**
     * Evaluates the spline at point t using the computed polynomial coefficients.
     * 
     * POLYNOMIAL EVALUATION:
     * For interval [x[i], x[i+1]], the cubic polynomial is:
     * S_i(t) = a[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])² + d[i]*(t-x[i])³
     * 
     * where splineCoeffs[i] = [d, c, b, a]:
     * - a[i] is the constant term (function value at knot x[i])
     * - b[i] is the coefficient of the linear term
     * - c[i] is the coefficient of the quadratic term  
     * - d[i] is the coefficient of the cubic term
     * 
     * @param t evaluation point
     * @return interpolated value
     */
    private double value(double t) {

        // TODO: I'm not sure if we should do extrapolation since values outside the knots are not meaningful/shouldn't exist in the tree.
        if (t < x[0]) {
            // Extrapolate linearly before first knot using first derivative
            double dx = t - x[0];
            return splineCoeffs[0][3] + splineCoeffs[0][2] * dx;
        }
        if (t > x[x.length-1]) {
            // Extrapolate linearly after last knot using first derivative
            int n = x.length - 1;
            double dx = t - x[n];
            return splineCoeffs[n-1][3] + splineCoeffs[n-1][2] * dx;
        }
        
        // Find the appropriate interval [x[i], x[i+1]] containing t, if outside the knots, use the first or last interval
        int i = findInterval(t);
        if (i < 0) {
            i = 0;
        }
        if (i >= x.length - 1) {
            i = x.length - 2;
        }
        
        // Evaluate cubic polynomial: S_i(t) = a + b*dx + c*dx^2 + d*dx^3
        double dx = t - x[i];
        return splineCoeffs[i][3] + splineCoeffs[i][2] * dx + splineCoeffs[i][1] * dx * dx + splineCoeffs[i][0] * dx * dx * dx;
    }
    
    /**
     * Evaluates the first derivative of the spline at point t.
     * 
     * DERIVATIVE EVALUATION:
     * For interval [x[i], x[i+1]], the derivative of the cubic polynomial is:
     * S'_i(t) = b[i] + 2*c[i]*(t-x[i]) + 3*d[i]*(t-x[i])^2
     * 
     * This is obtained by differentiating:
     * S_i(t) = a[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])^2 + d[i]*(t-x[i])^3
     * 
     * @param t evaluation point
     * @return interpolated derivative value
     */
    private double derivative(double t) {
        if (t < x[0]) {
            // Extrapolate derivative before first knot (constant derivative)
            return splineCoeffs[0][2];
        }
        if (t > x[x.length-1]) {
            // Extrapolate derivative after last knot (constant derivative)
            return splineCoeffs[x.length-2][2];
        }
        
        // Find the appropriate interval [x[i], x[i+1]] containing t
        int i = findInterval(t);
        if (i < 0) {
            i = 0;
        }
        if (i >= x.length - 1) {
            i = x.length - 2;
        }
        
        // Evaluate derivative: S'_i(t) = b[i] + 2*c[i]*dx + 3*d[i]*dx²
        double dx = t - x[i];
        return splineCoeffs[i][2] + 2.0 * splineCoeffs[i][1] * dx + 3.0 * splineCoeffs[i][0] * dx * dx;
    }

    /**
     * Gets value at time t.
     * 
     * @param t time point
     * @return value at time t
     */
    public double getValue(double t) {
        return value(t);
    }

    /**
     * Gets value at time t using precomputed grid points for efficiency.
     * 
     * @param t time point
     * @return value at closest grid point
     */
    public double getValueAtGridPoint(double t) {
        // Clamp t to grid range
        if (t <= gridStart) {
            return gridValues[0];
        }
        if (t >= gridEnd) {
            return gridValues[gridValues.length - 1];
        }
        
        // Find closest grid point
        int index = (int) Math.round((t - gridStart) / gridStep);
        index = Math.max(0, Math.min(index, gridValues.length - 1));
        
        return gridValues[index];
    }

    /**
     * Gets derivative at time t.
     * 
     * @param t time point
     * @return derivative at time t
     */
    public double getDerivative(double t) {
        return derivative(t);
    }

    /**
     * Gets derivative at time t using precomputed grid points for efficiency.
     * Returns the derivative value at the closest grid point.
     * 
     * @param t time point
     * @return derivative at closest grid point
     */
    public double getDerivativeAtGridPoint(double t) {
        // Clamp t to grid range
        if (t <= gridStart) {
            return gridDerivatives[0];
        }
        if (t >= gridEnd) {
            return gridDerivatives[gridDerivatives.length - 1];
        }
        
        // Find closest grid point
        int index = (int) Math.round((t - gridStart) / gridStep);
        index = Math.max(0, Math.min(index, gridDerivatives.length - 1));
        
        return gridDerivatives[index];
    }


    
    /**
     * Finds the interval index i such that x[i] <= t < x[i+1].
     * 
     * @param t evaluation point
     * @return interval index, or -1 if t < x[0], or x.length-1 if t >= x[x.length-1]
     */
    private int findInterval(double t) {
        if (t < x[0]) return -1;
        if (t >= x[x.length-1]) return x.length - 1;
        
        // Binary search for efficiency
        int left = 0;
        int right = x.length - 1;
        
        while (left < right) {
            int mid = (left + right) / 2;
            if (x[mid] <= t && t < x[mid + 1]) {
                return mid;
            } else if (t < x[mid]) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        
        return left;
    }
    
    /**
     * Gets the number of knots.
     * 
     * @return number of knots
     */
    public int getKnotCount() {
        return x.length;
    }
    
    /**
     * Gets the knot at index i.
     * 
     * @param i index
     * @return knot value
     */
    public double getKnot(int i) {
        return x[i];
    }
    
    /**
     * Gets the function value at knot i.
     * 
     * @param i index
     * @return function value
     */
    public double getValue(int i) {
        return y[i];
    }
    
    /**
     * Gets the number of grid points.
     * 
     * @return number of grid points
     */
    public int getGridPointCount() {
        return gridPoints.length;
    }
    
    /**
     * Gets the grid step size.
     * 
     * @return grid step size
     */
    public double getGridStep() {
        return gridStep;
    }
    
    /**
     * Gets the grid start time.
     * 
     * @return grid start time
     */
    public double getGridStart() {
        return gridStart;
    }
    
    /**
     * Gets the grid end time.
     * 
     * @return grid end time
     */
    public double getGridEnd() {
        return gridEnd;
    }
}
