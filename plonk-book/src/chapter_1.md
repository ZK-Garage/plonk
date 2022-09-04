# Chapter 1: Background Material

## Groups
A binary operation $*$ on a set $G$ is a mapping from $G\times G$ to $G$, which associates to elements $x$ and $y$ of $G$ a third element $x*y$ of $G$. 

#### Definition: 
A group $(G, *)$ consists of a set $G$ together with a binary operation $*$ for which the following properties are satisfied: 

* Associativity: $(x*y)*z=x*(y*z),\; \forall x,y,z\in G $ 
* Neutral element: $\exists! \; e \in G,\; e*x=x=x*e,\; \forall x\in G$
* Inverse element: $\forall x\in G, \; \exists! \; x'\in G,\; x*x'=e=x'*x$ where $e$ is the neutral element of $G$. 

A group $G$ is Abelian (or commutative) if: $x*y=y*x,\forall x,y\in G$

### Cyclic groups:
 
#### Definition: 
A group $G$ is said to be cyclic, with generator $x$, if every element of $G$ is of
the form $x^n$ for some integer $n$.

## Fields
 
A field $F$ is a set with two binary operations $+$ and $*$ that satisfies the following field axioms:
 
* Closure under addition:  $\forall x,y\in F,\; x+y\in F$
* Closure under multiplication:  $\forall x,y\in F,\; x*y\in F$
* Additive inverses:  $\forall x\in F,  y\in F$ such that  $x+y=0$
* Multiplication inverses:  $\forall x\in F$ such that $x\neq 0, \; \exists \; y\in F$ such that $x*y=1$, 
$y$ is called the multiplicative inverse of $x$ and is denoted $x^{-1}$ or $\dfrac{1}{x}$
* The distributive law:  $\forall x,y,z\in F,\; x*(y+z)=x*y+x*z$

### Finite Fields:


A finite field is a field $F_q$ with a finite number of elements. The order of a finite field $$|F_q|=q=p^k$$ 
for some integer $k\geq 1$ and $p$ prime equals the number of elements in the field.

## Elliptic curves
For ARK-PLONK, we use pairing friendly elliptic curves defined over very large finite fields $|F_p|≈2^{256}$. In the following we will explain what is an elliptic curve, pairings and elliptic curves over finite fields.

### Definition
An Elliptic curve $E$ is a mathematical object defined over a field $F$ and generally expressed in the following Weierstrass form:
$$y^2=x^3+ax+b$$ 
for some $a,b\in F_q$ where $(x,y)$ are called "affine coordinates" 

### Elliptic curves over finite fields

We will mainly focus on elliptic curves over finite fields since they are the ones used for cryptographic applications. An elliptic curve over a finite field $F_q$ is an abelian group $G$ with a finite number of points $n$ such that $n=|G|$ (order of the group $G$). 

The group on which the elliptic curve is defined is the one used for many cryptographic protocols. An elliptic curve over a finite field is represented in the projective plane, such a curve will have an additional point at infinity $O$ which serves as the identity of the group. There are several ways to define the curve equation, but for the purpose of doing the group law arithmetics, let $y^2=x^3+b$ for some constant $b\in F_q$.

### Group law

We can define an abelian group $G$ over elliptic curves as follows:

* The elements of the group are the points of an elliptic curve
* The identity element is the point $O$ defined as $(0,1,0)$
* The inverse of a point $P=(x,y,z)$  is the point $-P=(x,-y,z)$
* Commutativity:  $P+Q=Q+P$
* Associativity: $P+(Q+R)=(P+Q)+R$


<figure>
<img src="images/images/EC.png"
     alt="EC arithmetic"
     style="float: left; margin-right: 10px;" />
<figcaption>Figure 1: Elliptic curve arithmetics</figcaption>
</figure>

We can identify four cases in the group law: point addition $P+Q=-R$
which basically draws the line that intersects both points and finds a third point $R$, the addition is simply $-R$ (negate the y-coordinate). The second case $P+Q+Q=O$ is basically when the addition of both points $P$ and $Q$ doesn’t result in a third point $R$ cause in this case the line passing through $P$ and $Q$ is tangent to the curve.  Let us assume that $Q$ is the tangency point, then $P+Q=-Q$. The third case $P+Q=O$ is when a point $P$ is being added to its negation $Q=-P$, their addition equals point at infinity $O$. The final case $P+P=O$ is when a point is being added to itself (called point doubling).

#### Point addition

We will first compute $P+Q$ in the affine form and then projective form. 

1. We have two points $P=(x_p,y_p)$ and $Q=(x_q,y_q)$ where $P\neq Q$. If $P=O$ then $P$ is the identity point which means $P+Q=Q$. Likewise if $Q=O$, then $P+Q=P$ \
Else, $P+Q=S$ for $S=(x_s,y_s)$ where $S=-R$ (the point represented in the graph) such that: $$x_s=\lambda^2-x_p-x_q$$ and $$y_s=\lambda(x_p-x_s)-y_p$$   where      $\lambda=\dfrac{y_q-y_p}{x_q-x_p}$
2. In the projective form, each elliptic curve point has 3 coordinates instead of 2 $(x,y,z)$ for $z\neq 0$ (for all points except the point at infinity). Using the projective form allows us to give coordinates to the point at infinity. It also speeds up some of the most used arithmetic operations in the curve. The forward mapping is given by $(x,y)\rightarrow (xz,yz,z)$ and reverse mapping is giving by $(x,y,z)(\dfrac{x}{z},\dfrac{y}{z})$. Let $P=(x_p,y_p,z_p)$ and $Q=(x_q,y_q,z_q)$ and $\dfrac{x_p}{z_p}\neq\dfrac{x_q}{z_q}$. 


<figure>
<img src="images/images/ECP.png"
     alt="projective"
     style="float: left; margin-right: 10px;" />
<figcaption>Figure 2: Graphical representation of the map with z =1</figcaption>
</figure>

By expanding the previous arithmetics, we have:

$\lambda=\dfrac{\dfrac{y_q}{z_q}-\dfrac{y_p}{z_p}}{\dfrac{x_q}{z_q}-\dfrac{x_p}{z_p}}=\dfrac{\dfrac{y_q}{z_q}-\dfrac{y_p}{z_p}}{\dfrac{x_q}{z_q}-\dfrac{x_p}{z_p}}*\dfrac{z_pz_q}{z_pz_q}$ 

$\lambda=\dfrac{y_qz_p-y_pz_q}{x_qz_p-x_pz_q}$ 


 $x_s=\lambda^2-\dfrac{x_p}{z_p}-\dfrac{x_q}{z_q}$ and     $y_s=\lambda(\dfrac{x_p}{z_p}-x_s)-\dfrac{y_p}{z_p}$

#### Point doubling

We will compute $P+P=2P$ in the affine form as follows:\
If $P=O$ then $2P=O$ \
Else  $P=(x,y)$:\
    &emsp;&nbsp;If $y=0$ then $2P=O$  
    &emsp;&nbsp;Else $2P=(x',y')$  

such that: the derivative $\lambda=\dfrac{dy}{dx}=\dfrac{3x^2}{2y},\; x'=\lambda^2-2x\;$   and $y'=\lambda(x-x')-y=\lambda^3+3\lambda x-y$

#### The discrete logarithm problem
The security of many cryptographic techniques depends on the intractability of the discrete logarithm problem.

#### Definition:
Let $G$ be a multiplicative group. The discrete logarithm problem (DLP) is:
Given $g,h\in G$ to find $a$, if it exists, such that $h=g^a$.


#### Elliptic Curve Discrete Logarithm Problem  (ECDLP)

Let $E$ be an elliptic curve of the Weierstrass form defined over a finite field $F_q$. Let $S$ and $T$ be two points in $E(F_q)$. Find an integer $m$ such that:
$$T=mS$$

The fastest method to solve the ECDLP problem in $E(F_q)$ is the [Pollard Rho](https://www.ams.org/mcom/1978-32-143/S0025-5718-1978-0491431-9/S0025-5718-1978-0491431-9.pdf) method which has exponential complexity $O(\sqrt{|G|})$. In order for this algorithm to be exponential, we need to define elliptic curves over very large fields $|F_p|≈2^{256}$ which is currently the case for ARK-PLONK.

### Pairings
Pairing based-cryptography is used in many cryptographic applications like signature schemes, key agreement, zero knowledge...etc.\
For example, pairings are used to create efficient circuit-based zero knowledge proofs.

#### Definition
A pairing $e$ is a bilinear map $〈·,·〉$  defined as: 
$$e:G_1\times G_2\rightarrow G_T$$
Such that $G_1,G_2$ and $G_T$ are abelian groups. The bilinear property means that:
$e(P+P',Q)=e(P,Q)+e(P',Q)$\
$e(P,Q+Q')=e(P,Q)+e(P,Q')$

​​From which, we deduce the following transformations:


$$e([a]P,[b]Q)=e(P,[b]Q)^a=e([a]P,Q)^b=e(P,Q)^{ab}=e([b]P,[a]Q)$$

for $P,P'\in G_1$ and $Q,Q'\in G_2$ and $a,b\in \Z$

### Pairing friendly curves

In zk-SNARK schemes such as PLONK, we need to manipulate very large polynomials in order to perform efficient multi-point evaluation/interpolation with fast fourier transforms. We therefore target a subfamily of curves that have:
* Optimal extension field towers.
* Simple twisting isomorphisms.
 
#### Montgomery curve

Montgomery curves were first introduced by Peter L. as a way to accelerate elliptic curve methods of factorization  and then became central to elliptic curve cryptography. Later, these same qualities also led to many efficient implementations of elliptic curve cryptosystems, most notably Bernstein’s Curve25519 software.


##### Definition: 
A Montgomery curve over $F_q$ is an elliptic curve defined as
$$E_{(A,B)}: By^2=x(x^2+Ax+1)$$
where $A$ and $B$ are parameters in $F_q$ satisfying $B\neq 0$ and $A^2\neq 4$. 
 
#### Edwards Curves and Twisted Edwards Curves

Edwards curves are a family of elliptic curves that uses Edwards form instead of the more well known Weierstrass form used by elliptic curves. Edwards curves provide better efficiency and security than a general elliptic curve. For example the [edwards25519](http://cr.yp.to/ecdh.html) curve is faster than [secp256k1](https://en.bitcoin.it/wiki/Secp256k1) curve and considered to be more secure according to the [safe curve security assessments](https://safecurves.cr.yp.to/twist.html). The edwards25519 is considered to be secure against [twist attacks](https://safecurves.cr.yp.to/twist.html) (small group and invalid group attacks) whereas the secp256k1 is considered to be vulnerable to those.

In ARK-PLONK, the user is free to choose which pairing friendly and embedded curve they wish to work with. Currently, we have only tested this functionality with the pairing friendly [Bls12-381](https://z.cash/blog/cultivating-sapling-new-crypto-foundations/) and [Bls12-377](https://docs.rs/ark-bls12-377/latest/ark_bls12_377/) which have embedded [jubjub](https://z.cash/technology/jubjub/) and [this twisted edwards](https://docs.rs/ark-ed-on-bls12-377/0.3.0/ark_ed_on_bls12_377/) curve but theoretically any curve that satisfies the criteria for the PLONK protocol should work fine (for example [BN254](https://docs.rs/ark-bn254/latest/ark_bn254/) and [babyjubjub](https://iden3-docs.readthedocs.io/en/latest/iden3_repos/research/publications/zkproof-standards-workshop-2/baby-jubjub/baby-jubjub.html) which are already implemented in arkworks and thus they should work fine).

##### Definitions:

#### Edwards curves:
An Edwards curve over $k$ is a curve: 

$$E: x^2+y^2=1+dx^2y^2-x^2$$ 
where $d\in k-\{0,1\}$ and $k$ is a field with $char(k)\neq 2$.

#### Twisted Edward Curves:

Fix a field $k$ with $char(k)\neq 2$ and distinct non-zero elements $a,d\in k$. The Twisted Edwards curve with coefficients $a$ and $d$ is the curve $$E_{E_{a,b}}:ax^2+y^2=1+dx^2y^2$$
An Edwards curve is a Twisted Edwards Curve with $a=1$.


#### Correspondence with Montgomery curves

Every Twisted Edwards curve is birationally equivalent to an elliptic curve in Montgomery form, and vice versa. 

#### Curve25519
Curve25519 is a Montgomery curve providing 128 bits of security defined as:    $$y^2=x^3+ax^2+x$$
over prime field $p$ where: $b=1$.\
The curve is birationally equivalent to a twisted Edwards curve used in the Ed25519 signature scheme.

#### Barreto-Naehrig curves(BN curves)
The [BN-curve](https://eprint.iacr.org/2005/133.pdf) is a pairing friendly elliptic curve built over a field $F_q$ for $q\geq 5$ that achieves both high security and efficiency and has optimal ate pairing. BN curves with 256-bit (for example BN256) were believed to provide a 128-bit security level, but due to recent research [exTNFS](https://link.springer.com/article/10.1007/s00145-018-9280-5) this number dropped to 100-bits of security.

#### Barreto-Lynn-Scott curves (BLS curves)

A BLS curve is a pairing over BLS curves that constructs optimal Ate pairings. BLS12-381 is optimal for zk-SNARKs at the 128-bit security level and is implemented by the zcash team. Bls12-381 has an embedded Jubjub curve. 
  
BLS-12 curves are a more efficient choice than BN curves in terms of optimal Ate pairings, and they have a better security level (Bls12-381 provides 128-bits security whereas BN256 provides only 100-bits security) 

**Jubjub curve:** is a twisted Edwards curve of the form $-x^2+y^2=dx^2y^2$ built over the BLS12-381 scalar field.


## Lagrange basis and polynomial interpolation

Polynomial interpolation is a process where a given set of points $(x_i, y_i)$, $i \in [n]$ allows us
to construct a polynomial $f(x)$ that passes through all of them. We will assume
that $x_i \neq x_j$ for all distinct $i,j$ pairs, otherwise there is a repeated pair or it is not
possible to construct the polynomial as it would have to take two different values at the same $x$-point.

Notice that for a set of 2 points, we can find a line that crosses both of them, for a
set of 3 points, a parabola, and in general, for a set of $n$ points there is a polynomial of degree
$n-1$ that contains all of them.


The Lagrange interpolation consists of 2 steps:

 1. Construct a **Lagrange basis**. This is a set of $n$ polynomials of degree $n-1$ that
    take the value 0 at all points of the set except one, where their value is 1. Expressed in a formula:
 $$
 L_i(x) =
 \begin{cases}
   0, & \text{if }\ x = {x_j}, j \in [n], j\neq i\\
   1, & \text{if }\ x = x_i
 \end{cases}
 $$
 The polynomials $L_i$ can be constructed in the following way:
 $$
 L_i(x) = \prod_{0 \leq j < n\text{,  }j\neq i}
 \frac{x-x_j}{x_i - x_j}
 $$
 Notice that this product has $n-1$ terms, and therefore results in a degree $n-1$ polynomial.

 2. Scale and sum the polynomials of the basis.
 $$
  f(x) = \sum_{i=0}^{n-1} y_i \cdot L_i(x)
 $$
 The properties of the Lagrange basis now allow us to scale each polynomial to its target vale --
 multiplying by $y_i$ and then add up all the terms.


The important observation we can extract from the Lagrange interpolation is that given a fixed set
of points $x_1,\dots,x_n$ (an evaluation domain) we can represent any polynomial of degree $d<n$
by its evaluations $f(x_i)$ at $d+1$ points in the set. As it turns out, this representation is
much more convenient than the usual coefficient representation as it provides a very simple and fast
way of computing sums and multiplication of polynomials.
However, the coefficient form is still useful for evaluating the polynomial at points outside the evaluation domain.

Switching between these two forms of representation is very useful. The coefficient form is
preferred when the polynomial must be evaluated at a random point (outside of the evaluation
domain). The evaluation form is better suited for operations between polynomials such as
addition, multiplication and exact quotients. The algorithm that allows us to efficiently switch
between representations is the Fast Fourier Transform (FFT).  This is an efficient algorithm for the
more general discrete Fourier Transform (DFT). It has a complexity of $\mathcal{O}(n \cdot log(n))$
with $n$ being the degree of the polynomial.

#### Fast Fourier Transform algorithm

The FFT is generally defined over the complex numbers but in the crypto context it is always used over
a finite field $\mathbb{F}$.  The only requisite for $\mathbb{F}$ is that it have a large
multiplicative subgroup $H$ of order $n=2^k$ for some $k \in \mathbb{N}$.  This subgroup $H$ will be
the evaluation domain and it will consist of the $n^{th}$ roots of unity
$$
H = \{ \omega, \omega^2, \dots, \omega^n \} =
\{ x \in \mathbb{F} | x^n -1 =0 \}
$$

## Commitment schemes
A commitment scheme $C$ is a protocol between two parties: a prover $P$ and a verifier $V$. The goal of such a scheme is to satisfy the following security properties:

* **Hiding:** $P$ should be able to commit to a value $m$ by encoding it using a key $P_K$ without $V$ learning any information about $m$.
* **Binding:** $P$ cannot “cheat” by sending a different key $P_K'$ which opens $C$ to a different value $m'$. 
 
This is very useful in various cryptographic applications including zero-knowledge proofs. A commitment scheme can either be interactive or non-interactive depending on the use case. Their security assumption also varies between perfect or computational security with respect to the hiding and binding properties. 
 
There are different combinations of these properties but the most famous ones are perfectly binding and computationally hiding commitment schemes, and computationally binding and perfectly hiding commitment schemes.\
In the first scheme, $P$ generates the public key and sends it to $V$. The perfectly binding means $P$ is unable to change the commitment value after it has been committed to. The computational hiding means the probability of $V$ being able to guess the commitment value is negligible.\
In the second scheme, $V$ generates the public key and sends it to $P$. The computational binding means the chance of being able to change the commitment is negligible. The perfectly hiding means that a commitment to a message $m$ reveals no information about $m$.
For a detailed explanation of those properties, check this [article](https://homepages.cwi.nl/~schaffne/courses/crypto/2014/papers/ComZK08.pdf). 
 
#### Definition:  
A non-interactive commitment scheme has the following three algorithms:
 
 1.  **Key generation (setup)** $(1^k)$: The algorithm key outputs a pair of keys $(P_K,V_K)$ that is sent to the prover and verifier respectively for a given security parameter $k$.
2.  **Commitment** $(P_K,m)$: The algorithm com takes as input the prover key $P_K$ and the message $m$ and outputs the commitment $C$ and an opening value $d$ which will be used by the verifier.


3.  **Verification** $(V_K,C,m,d)$: The algorithm takes the verification key $V_K, C, m$ and $d$ as input and outputs yes or no depending on whether the verification is successful.
 
### Kate Polynomial commitment scheme (PCS)

Plonk uses Polynomial commitments, in particular Kate commitments (KZG10) to construct a fully-succinct zk-SNARK protocol with a constant proof size and verification time. KZG10 allows the prover to commit to a set of polynomials and to then show correct evaluation of these polynomials in a set of points with only one small proof. 
KZG is not the only PC scheme that could be used to construct ZK Proof systems but is the best one in terms of succinctness compared to other schemes like FRI, DARK and IPA where the proof size increases with the increase of the circuit scale. For example, the DARK scheme has logarithmic proof size and linear verification time. In terms of security, the FRI-based ZK-STARKs algorithm provides both plausible quantum-resistance and does not require any trusted setups. KZG however uses elliptic curve pairings, which are not quantum resistant and requires a third-party trusted setup.

We will first provide a formal definition of the original version of Kate commitments in which the prover commits to a polynomial $f$ and then performs an evaluation opening at a single point. We will then introduce the batched version, which is the one actually used in Plonk. 

#### Terminology:

Let $\phi\in F_p[X]$ be a polynomial of degree $\leq t$ in $F_p[X]$ where $F_p[X]$ is a polynomial ring over the finite field $F_p$ of prime order $p$. Let $G$ be an elliptic curve group of prime order $p$ such that there exists a symmetric bilinear pairing $e:G\times G\rightarrow G_T$ into the multiplicative group $G_T$ and for which the t-Strong Diffie-Hellman (t-SDH) Assumption holds. Let $g$ be a generator of  $G$. We will be using the multiplicative notation of the DLP assumption explained earlier instead of the additive notation for the EC group. We will use the following notation: $[x]=g^x​​$



 
 
#### Definition: 
A polynomial  commitment  scheme  consists  of  six  algorithms:\
**Setup, Commit, Open, VerifyPoly, CreateWitness, and VerifyEval**.

1. **Setup:** This step is run by a trusted or distributed authority and generates a commitment pair of a public and secret key $(P_K,S_K)$ as follows: $$⟨\delta,[1], [\alpha], [\alpha^2],.........,[\alpha^t]⟩$$  where $\delta=(e,G,G_t)$ is a bilinear pairing group, $⟨[1], [\alpha], [\alpha^2],.........,[\alpha^t]⟩\in G^{t+1}$ and $S_K=\alpha$  such that $\alpha \in F_p$. 
Note that $S_K$ is not required in the rest of the scheme.


2. **Commit:** The prover computes a commitment $C$ which hides the polynomial $\phi(x)$ as follows: $$     C=[\phi(\alpha)] =\prod_{j=0}^{t}([\alpha^j])^{\phi_j}$$

Where $\phi(\alpha)$ is an evaluation of the polynomial $\phi$ at point $\alpha$, $\phi(x)=\sum_{j=0}^{t}\phi_jx^j$ and $t=deg(\phi)$. 
 
3. **Open:** The prover outputs the polynomial $\phi(x)$ used to create the commitment.


4. **VerifyPoly:** The verifier checks if $C$ is indeed a commitment to $\phi(x)$ as verifying if $C$ is equal to $[\phi(\alpha)]$. The verifier is able to do so since: $$[\phi(\alpha)]=\prod_{j=0}^{t}([\alpha^j])^{\phi_j}$$ and $⟨[1], [\alpha], [\alpha^2],.........,[\alpha^t]⟩$ is known to the verifier.


5. **CreateWitness:** The prover outputs $<\beta,\phi(\beta),w_{\beta}>$, where $w_{\beta}$ is a witness for the evaluation $\phi(\beta)$ of $\phi(x)$ at index $\beta\neq \alpha$ chosen by the verifier. $$w_{\beta}=g^{\Psi_{\beta}(\alpha)}$$  and $$\Psi_{\beta}(x)=\dfrac{\phi(x)-\phi(\beta)}{x-\beta}(x)$$

The construction of $\Psi$ is based on an algebraic property of polynomials which guarantees that $\psi(x)-\psi(\beta)$ is divisible by $x- \beta$ for any $\beta\in F_p$.
 
6. **VerifyEval:** In this last step, the verifier checks if $\phi(\beta)$ is indeed the evaluation of $\phi(x)$ at index $\beta$ by checking if the following holds: $$e(C,g)=e(g^{\Psi_{\beta}(\alpha)},g^{\alpha}.g^{-\beta})e(g,g)^{\phi(\beta)}$$
   

## Zero Knowledge Proof systems

Before we start talking about zero knowledge proof systems, let's first define what a proof system is:\
A proof system is a protocol by which one party (prover) wants to convince another party (verifier) that a given statement is true. 

### Zero knowledge proof system:
In zero-knowledge proofs, the prover convinces the verifier about the truthfulness of the statement without revealing any information about the statement itself.
 
#### Properties
A zero-knowledge proof needs to fulfill each of the following properties to be fully described:
* **Completeness:**  An honest prover is always able to convince the verifier of the truthfulness of their claim.
* **Soundness:** If the prover’s claim is false (malicious prover), the verifier is not convinced.
* **Zero knowledge:** The proof should not reveal any information to the verifier beyond the truthfulness of the given claim.

#### Types of Zero knowledge proofs

There are two types of zero knowledge proofs (interactive and non interactive ones)
Interactive zero-knowledge proof systems were first introduced in 1985 by Goldwasser, Micali and Rackoff. Non-interactive schemes were introduced later on by Blum et al. The main difference between both schemes is that interactive proofs require interaction between both parties which means both have to be online in order to do so; this can be seen as inconvenient, especially for modern cryptography applications, while non-interactive proofs need a shared setup preprocessing phase instead. The shared setup phase will allow the participating parties to know which statement is being proved and what protocol is being used.


#### Interactive Zero Knowledge Proofs 

A prover $P$ has a secret $s$ and correctly responds to challenges to convince a verifier $V$ it has knowledge of $s$ using rounds of interaction between the two parties.
 
#### Non interactive zero knowledge proofs (NIZK)

Non-interactive zero-knowledge proofs, also known as NIZKs are another type of zero-knowledge proof which require no interaction between the prover and the verifier. 
In order to transform interactive proofs into NIZK proofs, cryptographers used the Fiat-Shamir heuristic hash function. This hash function allows one to compute the verifier challenges and offer very efficient NIZK arguments that are secure in the random oracle model. More recent works have started using bilinear groups to improve efficiency.


The two major types of NIZK proofs are zkSNARKs and zkSTARKs; zkSNARKs are based on elliptic curve cryptography and need a trusted setup phase, whereas zkSTARKs rely on hash functions and do not have any trusted setup. 
We will focus on zkSNARKs since PLONK is in this category.

zkSNARKs stands for zero-knowledge Succinct Non-interactive ARguments of Knowledge. 
  
- Succinct: proof length needs to be short
- Non-interactive: needs to be verifiable in a short amount of time
- ARKs: need to show that we know an input (witness) which yields to a certain computation. 

zkSNARKs cannot be applied to any computational problem directly; rather, you have to convert the problem into the right “form” for the problem to operate on. 

 

 



