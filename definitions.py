def adjoint(M):
    return M.transpose().conjugate();

# The complex conjugate of a transposed matrix is called its Hermitian conjugate, denoted by a dagger †.
def hermitianConjugate(M):
    return adjoint(M);

# Get a ⟨k| from a given |k⟩ or v.v.
def dual(state):
    return adjoint(state);

# Normalize a state
def normalize(state):
    return state/sqrt((state.conjugate().transpose() * state)[0][0]);


ket0 = matrix([[1], [0]]);             # |0⟩
bra0 = dual(ket0);                     # ⟨0|

ket1 = matrix([[0], [1]]);             # |1⟩
bra1 = dual(ket1);                     # ⟨1|

ketPlus = 1/sqrt(2) * (ket0 + ket1);   # |+⟩
braPlus = dual(ketPlus);               # ⟨+|

ketMinus = 1/sqrt(2) * (ket0 - ket1);  # |-⟩
braMinus = dual(ketMinus);             # ⟨-|

ketI = 1/sqrt(2) * (ket0 + i * ket1);  # |i⟩
braI = dual(ketI);                     # ⟨i|

ketINeg = 1/sqrt(2) * (ket0 - i * ket1); # |-i⟩
braINeg = dual(ketI);                    # ⟨-i|

vectorStates = [ket0, ket1, ketPlus, ketMinus, ketI, ketINeg];

ket00 = ket0.tensor_product(ket0); # |00⟩
ket01 = ket0.tensor_product(ket1); # |01⟩
ket10 = ket1.tensor_product(ket0); # |10⟩
ket11 = ket1.tensor_product(ket1); # |11⟩

ket000 = ket00.tensor_product(ket0); # |000⟩
ket001 = ket00.tensor_product(ket1); # |001⟩
ket010 = ket01.tensor_product(ket0); # |010⟩
ket011 = ket01.tensor_product(ket1); # |011⟩
ket100 = ket10.tensor_product(ket0); # |100⟩
ket101 = ket10.tensor_product(ket1); # |101⟩
ket110 = ket11.tensor_product(ket0); # |110⟩
ket111 = ket11.tensor_product(ket1); # |111⟩

# Identity
I = matrix([[1, 0], [0, 1]])
for state in vectorStates:
    assert(I * state == state)

# Pauli X operator flips about the x-axis and represents negation.
X = (ket0 * bra1) + ket1*bra0;
NOT = X;
SigmaX = X;
σ_x = X;
assert(X * ket0 == ket1)
assert(X * ket1 == ket0)

# Pauli Y operator (or $\sigma_y$) flips about the y-axis.
Y = matrix([[0, -i], [i, 0]]);
SigmaY = Y;
σ_y = Y;
assert(Y * ket0 == i * ket1)
assert(Y * ket1 == -i * ket0)

# Pauli Z operator (or $\sigma_y$) flips the state about the z-axis.
Z = matrix([[1, 0], [0, -1]]);
SigmaZ = Z;
σ_z = Z;
assert(Z * ket0 == ket0)
assert(Z * ket1 == - ket1)

# S operator rotates the state about the z-axis by 90°
S = matrix([[1, 0], [0, i]]);
assert(S**2 == Z)
assert(S * ket0 == ket0);
assert(S * ket1 == i * ket1);

T = matrix([[1, 0], [0, exp(i * pi / 4)]]);
assert(S == T**2);
assert(T * ket0 == ket0);
assert(T * ket1 == exp(i * pi/4) * ket1);

def R(r):
    return matrix([[1, 0], [0, exp(i * r)]])

assert(R(𝜋) == Z)
assert(R(𝜋/2) == S)
assert(R(𝜋/4) == T)

# Hadamard operator
H = 1/sqrt(2) * matrix([[1, 1], [1, -1]]);

# The Hadamard operator puts a state into the superposition of states
assert(H * ket0 == (1/sqrt(2) * (ket0 + ket1)))

assert(H * X * H == Z)
assert(H * Z * H == X)
assert(H * Y * H == -Y)

assert(hermitianConjugate(H) == H)
assert(H**2 == I)

# An opertor that satisfied $U{\dagger}U = I$ is called unitary.
assert(hermitianConjugate(H) * H == I)

SWAP = matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]);
assert(SWAP * ket01 == ket10)

# The first quibit is the control quibit and the second is the target. If
# control is |0⟩ then we do nothing. Otherwise we apply NOT.
CNOT = matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]);
assert(CNOT * ket10 == ket11)

CZ = matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]);

# Both control qubits have to be |1⟩ for the target to be manipulated.
CCNOT = matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 1, 0],
                ]);
Toffoli = CCNOT;
assert(CCNOT * ket110 == ket111);

# If the control is |1⟩ we swap the remaining two.
CSWAP = matrix([[1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
               ]);
Fredkin = CCNOT;
assert(CSWAP * ket110 == ket101)

from sage.plot.plot3d.plot3d import axes
from sage.plot.plot3d.shapes import Text

# Convert spherical coordinates to cartesian coordinates
def spherical_to_cartesian_coordinates(point):
    theta = point[0];
    phi = point[1];
    return (
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    );

# Convert a state to spherical coordinates
def state_to_spherical_coordinates(state):
    a = state[0][0];
    b = state[1][0];
    
    # The state must be normalized
    assert(a*conjugate(a) + b*conjugate(b) == 1)
    # assert(sum(list(map(lambda c: c*conjugate(c), ketI.coefficients()))) == 1)
    
    # Remove the global phase
    a = abs(a);
    b = b / e^(i*(arg(a)));
    
    forget();
    phi, theta = var('phi, theta');
    assume(phi, 'real');
    # assume(phi >= 0);
    # assume(phi < 2 * pi);
    assume(theta, 'real');
    assume(theta >= 0);
    assume(theta <= pi);
    solutions_theta = solve(cos(theta/2) == a, theta, solution_dict=True);
    if(len(solutions_theta) < 1):
        raise Exception('No solutions found for real valued 0 ≤ Θ < 2π');
    solved_theta = solutions_theta[0][theta];

    solutions_phi = solve(
        e^(i * phi) * sin(solved_theta/2) == b,
        phi,
        algorithm='maxima',
        solution_dict=True
    );
    if(len(solutions_phi) < 1):
        solutions_phi = solve(
            e^(i * phi) * sin(solved_theta/2) == b,
            phi,
            algorithm='sympy',
            domain='real',
            solution_dict=True
        );
        if(len(solutions_phi) < 1):
            raise Exception('No real solutions found for ɸ');
    solved_phi = solutions_phi[0][phi];

    return (solved_theta, solved_phi);

# Convert a state to cartesian coordinates
def state_to_cartesian_coordinates(state):
    return spherical_to_cartesian_coordinates(
        state_to_spherical_coordinates(state)
    );

# Get a vector for display on a Bloch sphere for a given state
def bloch_vector(state, rgbcolor=(1, 0, 0), label=''):
    state_cartesian = state_to_cartesian_coordinates(state)
    state_cartesian_numeric = (
        n(state_cartesian[0]),
        n(state_cartesian[1]),
        n(state_cartesian[2]),
    );
    vec_state = arrow3d((0, 0, 0), state_cartesian_numeric, radius=0.01, alpha=1.0, rgbcolor=rgbcolor)
    
    if label != '':
        text_label = Text(
            label,
            rgbcolor=(0, 0, 0)
        ).translate(
            state_cartesian_numeric[0] + 0.1,
            state_cartesian_numeric[1] + 0.1,
            state_cartesian_numeric[2] + 0.1
        );
        return vec_state + text_label;
    
    return vec_state;

# Get a 3d Graphic of a Bloch sphere representing the given state
def bloch(state, vector_rgbcolor=(1, 0, 0)):
    a = axes(1, radius=0.1, color='black', alpha=0.5);

    circle_real = parametric_plot3d((sin, cos, 0), (0, 2*pi), rgbcolor=(0, 0, 0), alpha=0.5);
    circle_imag = parametric_plot3d((sin, 0, cos), (0, 2*pi), rgbcolor=(0, 0, 0), alpha=0.5);

    vec_i = arrow3d((0, 0, 0), (0, -1, 0), radius=0.002, rgbcolor=(0, 0, 0), alpha=0.5);
    vec_1 = arrow3d((0, 0, 0), (0, 0, -1), radius=0.002, rgbcolor=(0, 0, 0), alpha=0.5);
    vec_neg = arrow3d((0, 0, 0), (-1, 0, 0), radius=0.002, rgbcolor=(0, 0, 0), alpha=0.5);

    label_0 = Text("|0⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(0, 0, 1 + 0.1);
    label_1 = Text("|1⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(0, 0, -(1 + 0.1));

    label_i = Text("|i⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(0, 1 + 0.1, 0);
    label_neg_i = Text("|-i⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(0, -(1 + 0.1), 0);

    label_plus = Text("|+⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(1 + 0.1, 0, 0);
    label_minus = Text("|-⟩", rgbcolor=(0, 0, 0), alpha=0.5).translate(-(1 + 0.1), 0, 0);

    return (
        circle_real + circle_imag
        + a
        + vec_i + vec_1 + vec_neg
        + label_0 + label_1 + label_i + label_neg_i + label_minus + label_plus
        + bloch_vector(state, rgbcolor=vector_rgbcolor)
    );