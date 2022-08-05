# PLONK

PLONK stands for Permutations over Lagrange-bases for Oecumenical Non Interactive Arguments of Knowledge. From hereon we shall simply write plonk.

Plonk is a general-purpose zero-knowledge proof scheme which solves a huge issue inherited in traditional zksnarks proof systems like Groth16: the _non-universal_ one-time trusted setup.
The trusted setup is a preprocessing phase which creates a Structured Reference String (SRS) that is available to both prover and verifier. The reason why an SRS is created is to prevent the prover from cheating and creating fake proofs and thus fulfilling the soundness property. The problem with a non-universal one-time trusted setup is that it’s a one-time event which means _every_ circuit needs a new SRS to be generated, which results in slow verification time. Plonk solves this problem by creating a single SRS that’s used for an unlimited number of arbitrary circuits (of a certain maximum size). This string is also updatable, which improves security.

## PLONK components

Plonk can be constructed as follow:

Program (code)→ Arithmetic circuit → Constraint systems → Permutation checks → Polynomial commitments (A batched version of KZG10)

We will use the following example to be able to explain each step and the transition from one to the other. For a more detailed explanation, please check the [original PLONK paper](https://eprint.iacr.org/2019/953.pdf) as it contains formal definitions and some interesting insights regarding efficiency etc.

**Problem definition:**

Prover wants to prove to Verifier that she knows the solution to the equation:
$$x^3+x+5=0$$

The goal is for Prover to evaluate the above function without revealing anything about the secret value $x$ (solution to the equation).


Prover creates a program to represent the problem in a function code, which then will be translated into an arithmetic circuit.


### Arithmetic circuit
This step transforms a program into an arithmetic circuit where two basic components are being used: wires and gates. Plonk uses fan-in two gates; therefore each gate has a left input, a right input, and an output. A circuit with $n$ gates will have $3n$ wires.

Plonk is gate-based instead of R1CS-based like some proof systems, such as Groth16. A primary difference between the two systems is in how they handle addition gates; in R1CS, addition gates are cheap since wires that go from an addition to a multiplication gate are not labeled, which is not the case for a gate-based system. The reason why plonk uses a gate-based system rather than an R1CS system is that the linear constraints (which are just wiring constraints) can be reduced to a permutation check. To better understand the advantages and disadvantages of each of those designs check this [article](https://hackmd.io/@aztec-network/plonk-arithmetiization-air#How-does-all-this-relate-to-R1CS).


The following circuit translates the previous equation $x^3+x+5=0$, giving 2 multiplication gates and 2 addition gates.

<p align="center">
  <img  alt="circuit" width="250"src="images/images/circuit.png" />

</p>



### Constraint system
The circuit is converted into a system of equations where the variables are the values on each of the wires, and there is one equation per gate. Let’s take our previous example:

$$\begin{equation*}
\begin{dcases}
\begin{align}
  a_1*b_1-c_1 &=0 & \\
  a_2*b_2-c_2 &=0 &\\
  a_3+b_3-c_3&=0 &\\
  a_4+b_4-c_4&=0
\end{align}
\end{dcases}
\end{equation*}$$

The final result is  $x^3+x+5-c_4=0$, which represents the program we wanted to solve for $c_4=x^3+x+5=0$.
The setup for each of those equations is of the form:
$$(Q_{L_i})a_i + (Q_{R_i})b_i +(Q_{O_i})c_i + (Q_{M_i})a_ib_i  + Q_{C_{i}} =0$$

for  $L$= left, $R$ = right, $O$= output, $M$= multiplication, $C$= constant
The arithmetic gates are modeled with the selector vectors:  $(Q_L,Q_R,Q_O,Q_C,Q_M)$

Each $Q$ value is a constant. We define $Q$ for an additive gate and a multiplicative one and a constant gate as follow:\
For an addition gate, we set:
$$Q_{L_i}=1, Q_{R_i}=1,Q_{M_i}=0,Q_{O_i}=-1,Q_{C_i}=0$$
For a multiplication gate, we set:
$$Q_{L_i}=0, Q_{R_i}=0,Q_{M_i}=1,Q_{O_i}=-1,Q_{C_i}=0$$

For a constant gate setting $a_i$ to some constant $x$, we set:
$$Q_{L_i}=1, Q_{R_i}=0,Q_{M_i}=0,Q_{O_i}=0,Q_{C_i}=-x$$


There are two types of constraints:


1. **Gate constraints:** Which represent the equations between the wires attached to the same gate. For example the equation $(1)$:  $$a_1*b_1-c_1=0 $$.


2. **Copy constraints:**  Plonk enforces copy constraints; these associate wires from the whole circuit that have the same value, for example the output of one gate would be associated with the input of its destination gate. These constraints are checked with a permutation argument. For example, $c_1$ is the output of equation $(1)$ and is also an input for another gate, so we copy that constraint into a new one $a_2$ and we claim equality $c_1=a_2$.


### Permutation checks

We introduce a permutation argument  used to assure the correct execution of the circuit. It allows to check the connection between different wires inside of the circuit and make sure that the output of a certain circuit is equal to the input of another for example $(c_2=b_3)$  where $c_2$ is the output of circuit $2$ and $b_3$ is the right input of circuit $3$.

Let $\phi_1,..........,\phi_k \in F_{<d}[X]$ and $\sigma :[kn]\rightarrow [kn] $ for $k$ = number of wires. We say for a set of polynomials $(g_1,.......,g_k) \in (F_{<d}[X])^k$ that $$(g_1,g_2,........,g_k)=\sigma(\phi_1,...........,\phi_k)$$ if the following holds:

$g_{(l)}=\phi_{(\sigma(l))}$  for $l \in [kn]$  where $\phi_{((j-1).n+i)}=\phi_j(w^i), g_{((j-1).n+i)}=g_j(w^i)$

for each $j\in [k], i\in [n]$

The prover will be able to select an appropriate $\sigma$ for a set of wire connections and the
verifier will be convinced that the connections are correct by checking the permutation argument on the set of polynomials with itself:

$$
 (\phi_1,\dots,\phi_k) = \sigma((\phi_1,\dots,\phi_k))
$$

Note that if this equality holds, we can substitute $(\phi_1,\dots,\phi_k)$ in the right hand side indefinitely as so:

$$
 (\phi_1,\dots, \phi_k)=\\
 \sigma((\phi_1,\dots, \phi_k)) =\\
 \sigma( \sigma((\phi_1, \dots, \phi_k)) ) = \dots
$$

and therefore the check effectively assures that

$$
 (\phi_1, \dots, \phi_k) = \sigma^i((\phi_1, \dots, \phi_3))
$$
for all $i \in \mathbb{N}$.




### KZG10 Batched commitments:

PLONK uses a batched kate commitment form in order to improve verifier efficiency by allowing for a parallel opening of commitments for each evaluation point possible.

Let’s take $t$ polynomials of degree $\leq d\;$. Let $F$ be a field of prime order. We denote by $F_{<d}[X]$ the set of polynomials over $F$ of degree $<d$. Let $G_1,G_2,G_t$ be groups of size $r$ and $e:G_1\times G_2\rightarrow G_t$ a bilinear pairing such that $e(g_1,g_2)=g_t$ with $g_1,g_2$ generators of $G_1$ and $G_2$ respectively.


##### Definition:
d-polynomial commitment scheme is a setting of $t$ polynomials $\phi_1,\phi_2,.......,\phi_t\in F_{<d}[x]$ of degree $d$, each such that $z_1,z_2,..........,z_t\in F$ are evaluation points for those polynomials. The alleged commitments to polynomials are $cm_1,cm_2,......,cm_t$ where $cm_i=com(\phi_i,srs)$ for $i\in[t]$ and alleged correct openings $s_1,s_2,.........,s_t$.


The commitment scheme has three steps, as follow:


* $gen(d)$: this step will generate a structured reference string $(srs)$ in a randomized way. The algorithm randomly chooses $x\in F$ and outputs

$$srs=([1]_1,[x]_1,[x^2]_1,[x^3]_1,.............[x^{d-1}]_1,[1]_2,[x]_2)$$
where   $[x]_1 =x.g_1$  and   $x_2=x.g_2$.


* $com(,srs)$: the commitment is computed as follows,
$$com(,srs):=[\phi(x)]_1$$


* $open:$ we present two scenarios:

   1. All evaluation points are equal $z_1=z_2=........=z_t=z$ $$open(\{cm_i\},\{z_i\},\{s_i\})$$ for $i \in [t]$

       a. Verifier sends a random $\gamma\in F$

       b. Prover computes $$h(x)=\sum_{i=1}^{t}\gamma^{i-1}.\dfrac{\phi_i(x)-\phi_i(z)}{x-z}$$ and then uses $srs$ to compute the commitment $W$ and send it to verifier $$W=h[(x)]_1$$

   c.  Verifier computes the following:
    $$F=\sum_{i\in[t]}\gamma^{i-1}.cm_i \;and\;      v=[\sum_{i\in[t]}\gamma^{i-1}.s_i]_1$$
    and accepts iff $$e(F-v,[1]_2).e(-W,[x-z]_2)=1$$

2. Let $z,z'$ be two distinct evaluation points and $t_1,t_2$ be the number of polynomials
We will describe the protocol when there are two distinct points among $z_1,\dots, z_t$. Let $z, z'$
be the distinct evaluation points and $t_1,t_2$ then number of polynomials $\{f_i\}_{i\in [t_1]},
\{f_i'\}_{i\in [t_2]}$ evaluated in $z, z'$ respectively.

Note that these protocols are not zero-knowledge. The notion of zero-knowledge is not even well defined for polynomial commitments. At the end, when we present the full protocol we
will add blinders to add the zero-knowledge property to the Plonk protocol.

 3. $open( \{cm_i\}_{i \in [t_1]},\{cm_i'\}_{i \in [t_2]}, z, z', \{s_i\}_{i \in [t_1]} \{s_i'\}_{i \in [t_2]})$
     1. Verifier sends random challenges $\gamma, \gamma' \in \mathbb{F}$
     2. Prover computes polynomials
        $$
         h(X) := \sum_{i=1}^{t_1} \gamma^{i-1} \frac{f_i(X) - f_i(z)}{X-z}\\
         h'(X):= \sum_{i=1}^{t_2} \gamma'^{i-1} \frac{f'_i(X) - f'_i(z')}{X-z'}
        $$
        and sends commitments $W=[h(x)]_1, W'=[h'(x)]_1$.
     3. Verifier chooses random $r' \in \mathbb{F}$ and computes
        $$
         F:= \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot cm_i -
         \left[ \sum_{i=1}^{t_1} \gamma^{i-1} s_i \right]_1\right) +\\
         r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot cm'_i -
         \left[ \sum_{i=1}^{t_2} \gamma'^{i-1} s'_i \right]_1\right)
        $$
    4. Verifier accepts iff
        $$
        e(F +z \cdot W + r'z'\cdot W', [1]_2) =
        e(W + r'\cdot W', [x]_2)
        $$

Extending the left side of the check we get
        $$
        e(F +z \cdot W + r'z'\cdot W', [1]_2) =\\
        e(  \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot cm_i -
         \left[ \sum_{i=1}^{t_1} \gamma^{i-1} s_i \right]_1\right) +
         r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot cm'_i -
         \left[ \sum_{i=1}^{t_2} \gamma'^{i-1} s'_i \right]_1\right)\\ +
         z \cdot \sum_{i=1}^{t_1} \gamma^{i-1} \frac{[f_i(x) - f_i(z)]_1}{x-z} +\\
         r'z' \cdot \sum_{i=1}^{t_2} \gamma'^{i-1} \frac{[f'_i(x) - f'_i(z')]_1}{x-z'}, [1]_2) =\\
        e(  \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot [f_i(x)]_1 -
         \left[ \sum_{i=1}^{t_1} \gamma^{i-1} f_i(z) \right]_1\right) +
         r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot [f'_i(x)]_1 -
         \left[ \sum_{i=1}^{t_2} \gamma'^{i-1} f'_i(z) \right]_1\right)\\ +
         z \cdot \sum_{i=1}^{t_1} \gamma^{i-1} \frac{[f_i(x) - f_i(z)]_1}{x-z} +\\
         r'z' \cdot \sum_{i=1}^{t_2} \gamma'^{i-1} \frac{[f'_i(x) - f'_i(z')]_1}{x-z'}, [1]_2) =\\
        e(  \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot
        [f_i(x) - f_i(z)]_1 \right) +
        r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot
        [f'_i(x) - f'_i(z)]_1 \right)\\ +
        z \cdot \sum_{i=1}^{t_1} \gamma^{i-1}
        \frac{[f_i(x) - f_i(z)]_1}{x-z} +\\
        r'z' \cdot \sum_{i=1}^{t_2} \gamma'^{i-1}
        \frac{[f'_i(x) - f'_i(z')]_1}{x-z'}, [1]_2) =\\
        e(  \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot
        [f_i(x) - f_i(z)]_1 \right) (1 + \frac{z}{x-z}) +\\
        r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot
        [f'_i(x) - f'_i(z)]_1 \right) (1 + \frac{z'}{x-z'}),
        [1]_2) =\\
        e(  \left( \sum_{i=1}^{t_1} \gamma^{i-1} \cdot
        [f_i(x) - f_i(z)]_1 \right) (\frac{x}{x-z}) +\\
        r' \cdot \left( \sum_{i=1}^{t_2} \gamma'^{i-1} \cdot
        [f'_i(x) - f'_i(z)]_1 \right) (\frac{x}{x-z'}),
        [1]_2)
        $$

From the right side we get
        $$
        e(W + r' \cdot W', [x]_2) =\\
        e( \sum_{i=1}^{t_1} \gamma^{i-1}
        \frac{[f_i(x) - f_i(z)]_1}{x-z}\\ +
        r' \cdot \sum_{i=1}^{t_1} \gamma'^{i-1}
        \frac{f'_i(x) - f'_i(z')}{x-z'} , [x]_2)
        $$






## PLONK protocol



**Common preprocessed input:**

The prover only uses the $\mathbb{G}_1$ part of the srs,
therefore the $\mathbb{G}_2$ part is part of the verifier
preprocessed input

$$
 n,\\
 srs \rightarrow ([x]_1, [x^2]_1,\dots,[x^{n+2}]_1),\\
 (q_M, q_L, q_R, q_O, q_C)_{i=1}^n \rightarrow
 \begin{cases}
     q_M(X) = \sum_{i=1}^n q_{Mi} L_i(X)\\
     q_L(X)= \sum_{i=1}^n q_{Li} L_i(X),\\
     q_R(X)= \sum_{i=1}^n q_{Ri} L_i(X),\\
     q_O(X)= \sum_{i=1}^n q_{Oi} L_i(X),\\
     q_C(X)= \sum_{i=1}^n q_{Ci} L_i(X),
 \end{cases}\\
 \sigma(X) \rightarrow
 \begin{cases}
     S_{\sigma_1}(X) = \sum_{i=1}^n \sigma(i) L_i(X),\\
     S_{\sigma_2}(X) = \sum_{i=1}^n \sigma(n+i) L_i(X),\\
     S_{\sigma_3}(X) = \sum_{i=1}^n \sigma(2n+i) L_i(X),\\
 \end{cases}\\
$$

The public input values will be written as part of the set of wires.
Having $l$ wires as public inputs: $$l, (w_i)_{i \in [l]}$$



### Prover Algorithm

The prover input is the set of values assigned to each wire:

**Prover input:** $(w_i)_{i \in [3n]}$

**Round 1** -- Commit to wire values.

Generate random blinding scalars $(b_1,\dots,b_9) \in \mathbb{F}_p$
Compute wire polynomials
$$
    a(X) = (b_1 X + b_2)Z_H(X) + \sum_{i=1}^n w_i L_i(X)\\
    b(X) = (b_3 X + b_4)Z_H(X) + \sum_{i=1}^n w_{n+i} L_i(X)\\
    c(X) = (b_5 X + b_6)Z_H(X) + \sum_{i=1}^n w_{2n+i} L_i(X)
$$

Compute commitments:
$$
    [a]_1 := [a(x)]_1,
    [b]_1 := [b(x)]_1,
    [c]_1 := [c(x)]_1
$$

Output $([a]_1, [b]_1, [c]_1).$

**Round 2** -- Permutation polynomial

Compute challenges $\beta, \gamma \in \mathbb{F}_p$:
$$
\beta = hash(Transcript,0),\\
\gamma = hash(Transcript,1)
$$

Compute permutation polynomial $z(X)$:
$$
    z(X) = (b_7  X^2 +b_8  X +b_9 )Z_H(X) +\\
    L_1(X) +
    \sum_{i=1}^n \left(
            L_{i+1}(X) \prod_{j=1}^i
    \frac{  (\omega_j +\beta\omega^{j-1}+\gamma)
            (\omega_{n+j} + \beta k_1 \omega^{j-1} +\gamma)
            (\omega_{2n+j}+ \beta k_2 \omega^{j-1} +\gamma)}
         {  (\omega_j +\beta \sigma(j) +\gamma)
            (\omega_{n+j} + \beta \sigma(n+j) +\gamma)
            (\omega_{2n+j}+ \beta \sigma(2n+j) +\gamma)}
            \right)
$$

Compute $[z]_1 := [z(x)]_1$

Output $[z]_1$

**Round 3**

Compute quotient challenge $\alpha \in \mathbb{F}_p$
$$
    \alpha = hash(Transcript)
$$

Compute quotient polynomial $t(X)$:
$$
    t(X) :=\\
    \frac{1}{Z_H(X)} (a(X)b(X)q_M(X) +
                      a(X)q_L(X) +
                      b(X)q_R(X) +
                      c(X)q_O(X) +
                      PI(X) +
                      q_C(X)) +\\
    \frac{\alpha}{Z_H(X)} ((a(X) + \beta X + \gamma)
                           (b(X) + \beta k_1 X + \gamma)
                           (c(X) + \beta k_2 X + \gamma)
                           z(X)) -\\
    \frac{\alpha}{Z_H(X)} ((a(X) + \beta S_{\sigma_1}(X) + \gamma)
                           (b(X) + \beta S_{\sigma_2}(X) + \gamma)
                           (c(X) + \beta S_{\sigma_3}(X) + \gamma)
                           z(X \omega)) +\\
    (z(X) - 1)L_1(X)\frac{\alpha^2}{z_H(X)}
$$

Note that all terms of the polynomial are divided by $Z_H(X)$. This can be done because if all the constraints hold the then the polynomials are 0 in all elements of H, and therefore, divisible by $Z_H(X)$

Split $t(X)$ into $<n$ degree polynomials
$t_{lo}(X), t_{mid}(X), t_{hi}(X)$ so that:

$$
t(X) = t_{lo}(X) + X^n t_{mid}(X) + X^{2n} t_{hi}(X) +
$$

Compute $[t_{lo}]_1 := [t_{lo}(x)]_1$,
        $[t_{mid}]_1 := [t_{mid}(x)]_1$,
        $[t_{hi}]_1 = [t_{hi}(x)]_1$

Output $([t_{lo}]_1, [t_{mid}]_1, [t_{hi}]_1)$

**Round 4**

Compute evaluation challenge $\zeta \in \mathbb{F}_p$:
$$
    \zeta = hash(Transcript)
$$

Compute opening evaluations at $\zeta$:
$$
    \bar{a}:= a(\zeta),
    \bar{b}:= b(\zeta),
    \bar{c}:= c(\zeta),\\
    \bar{s}_{\sigma_1} := S_{\sigma_1}(\zeta),
    \bar{s}_{\sigma_2} := S_{\sigma_2}(\zeta),\\
    \bar{t} := t(\zeta),
    \bar{z}_{\omega} := z(\omega\zeta),
$$


Compute linearisation polynomial $r(X)$:
$$
r(X) =
(
    \bar{a}\bar{b} \cdot q_M(X) +
    \bar{a} \cdot q_L(X) +
    \bar{b} \cdot q_R(X) +
    \bar{c} \cdot q_O(X) +
     q_C(X)
    ) + \\
    (
    (\bar{a} + \beta \zeta + \gamma)
    (\bar{b} + \beta k_1 \zeta + \gamma)
    (\bar{c} + \beta k_2 \zeta + \gamma)
    ) \cdot z(X) \alpha - \\
    (
    (\bar{a} + \beta \bar{s}_{\sigma_1} + \gamma)
    (\bar{b} + \beta \bar{s}_{\sigma_2} + \gamma)
    \beta\bar{z}_\omega \cdot S_{\sigma_3}(X)
    )\alpha + \\
    z(X)L_1(\zeta)\alpha^2
$$

Compute linearisation evaluation at $\zeta$:
$$
    \bar r:= r(\zeta)
$$

Output ($\bar{a}, \bar{b},\bar{c},
        \bar{s}_{\sigma_1}, \bar{s}_{\sigma_2},
        \bar{z}_\omega, \bar{t}, \bar{r}$)

**Round 5**

Compute opening challenge $v \in \mathbb{F}_p$:
$$
    v = hash(Transcript)
$$

Compute opening proof polynomial $W_\zeta(X)$:

$$
    W_\zeta(X)=
    \frac{1}{X - \zeta}
    \left(
        (
        t_{lo}(X) +
        \zeta^n t_{mid}(\zeta) +
        \zeta^{2n} t_{hi}(\zeta)-
        \bar{t}
        ) +\\
        v(r(X) - \bar r) +\\
        v^2(a(X) - \bar a) +\\
        v^3(b(X) - \bar b) +\\
        v^4(c(X) - \bar c) +\\
        v^5(S_{\sigma_1} - \bar s_{\sigma_1}) +\\
        v^6(S_{\sigma_2} - \bar s_{\sigma_2})
    \right)
$$

Compute opening proof polynomial

$$
    W_{\zeta \omega}(X) =
    \frac{(z(X) - \bar z_\omega)}{X - \zeta\omega}
$$

Compute $[W_\zeta]_1:=[W_\zeta(x)],
[W_{\zeta \omega}]_1:=[W_{\zeta \omega}(x)]$


Return
$$
    \pi_{SNARK} =
    ([a]_1, [b]_1, [c]_1, [z]_1,
    [t_{lo}]_1, [t_{mid}]_1, [t_{hi}]_1,
    [W_\zeta]_1, [W_{\zeta \omega}]_1, \\
    \bar{a}, \bar{b}, \bar{c},
    \bar s_{\sigma_1}, \bar s_{\sigma_2},
    \bar{r}, \bar z_\omega
    )
$$

Compute multipoint evaluation challenge $u \in \mathbb{F}_p$:
$$
    u := hash(Transcript)
$$


### Verifier Algorithm


**Verifier preprocessed input:**

$$[q_M]_1 := [q_M(x)]_1,[q_L]_1 := [q_L(x)_1,[q_R]_1 := [q_R(x)]_1,$$
$$[q_O]_1 := [q_O(x)]_1,[q_C]_1 := [q_C(x)]_1,[s_{\sigma_{1}}]_1 := [S_{\sigma_{1}}(x)],$$
$$[s_{\sigma_{2}}]_1 := [S_{\sigma_{2}}(x)],[s_{\sigma_{3}}]_1 := [S_{\sigma_{3}}(x)],([1]_2, [x]_2)$$

**Verifier input:** $(w_i)_{[l]}, \pi_{SNARK}$

1. Validate $([a]_1, [b]_1, [c]_1, [z]_1,
              [t_{lo}]_1, [t_{mid}]_1, [t_{hi}]_1,
              [W_z]_1, [W_{z \omega}]_1 ) \in \mathbb{G}_1$

2. Validate $(\bar a, \bar b, \bar c,
            \bar s_{\sigma_1}, \bar s_{\sigma_2},
             \bar z, \bar z_\omega) \in \mathbb{F}_p^7$

3. Validate $(w_i)({i \in [l]}) \in \mathbb{F}_p^l$

4. Compute challenges from common input public input and the elements of $\pi_{SNARK}$: $\beta, \gamma, \alpha, \zeta, v, u \in \mathbb{F}_p$

5. Compute $Z_H(\zeta) = \zeta^n -1$

6. Compute $L_1(\zeta) = \frac{\omega (\zeta^n - 1)}{n (\zeta - \omega)}$

7. Compute Public input polynomial evaluation
 $PI(\zeta) = \sum_{i\in[l]} w_i L_i(\zeta)$

8. Compute quotient polynomial evaluation
$$
    \bar t =
    \frac{
        \bar r +
        PI(\zeta) -
        (
        (\bar a + \beta \bar s_{\sigma_1} + \gamma)
        (\bar b + \beta \bar s_{\sigma_2} + \gamma)
        (\bar c + \gamma) \bar z_\omega
        ) \alpha -
        L_1(\zeta) \alpha^2
    }{
        Z_H(\zeta)
    }
$$

9. First part of batched polynomial commitment
 $[D]_1 := v \cdot [r]_1 + u\cdot [z]_1$:
$$ v (\bar a \bar b [q_M]_1 +\bar a [q_L]_1 +\bar b [q_R]_1 +\bar c [q_O]_1 +[q_C]_1 )+$$
$$[z]_1 ((\bar a + \beta \zeta + \gamma)(\bar b + \beta k_1 \zeta + \gamma)(\bar c +
\beta k_2 \zeta + \gamma)\alpha v + L_1(\zeta)\alpha^2 v +u)$$
$$ - (\bar a + \beta \bar s{\sigma_1} \zeta + \gamma)(\bar b + \beta \bar s_{\sigma_2} \zeta + \gamma)\alpha v\beta \bar z_\omega[s_{\sigma_3}]_1$$

10. Compute the fully batched polynomial commitment $[F]_1$:
$$
    [F]_1 :=
        [t_{lo}]_1 +
        \zeta^n [t_{mid}]_1 +
        \zeta^{2n} [t_{hi}]_1  +
        [D]_1 +
        v^2 \cdot [a]_1 +
        v^3 \cdot [b]_1 +
        v^4 \cdot [c]_1 +
        v^5 \cdot [s_{\sigma_1}]_1 +
        v^6 \cdot [s_{\sigma_2}]_1
$$

11. Compute group-encoded bath evaluation $[E]_1$:
$$
    [E]_1 :=
    \left[
        \bar t +
        v \bar r +
        v^2 \bar a +
        v^3 \bar b +
        v^4 \bar c +
        v^5 \bar s_{\sigma_1} +
        v^6 \bar s_{\sigma_2} +
        u \bar z_\omega
    \right]_1
$$

12. Batch validate all evaluations
$$
    e(
        [W_\zeta]_1) +
        u \cdot [W_{\zeta \omega}]_1),
        [x]_2
    )
    \stackrel{?}{=}
    e(
        \zeta \cdot [W_\zeta]_1) +
        u \zeta \omega \cdot [W_{\zeta \omega}]_1) +
        [F]_1 - [E]_1,
        [1]_2
    )
$$

