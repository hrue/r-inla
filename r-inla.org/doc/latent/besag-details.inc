%%\documentclass{article}
%%\usepackage{amsmath}
%%\begin{document}

\section*{Details on the implementation}

This gives some details of the implementation, which depends on the
following variables
\begin{description}
\item[nc1] Number of connected components in the graph with size 1.
    These nodes, \emph{singletons}, have no neighbours.
\item[nc2] Number of connected components in the graph with size
    $\ge2$.
\item[scale.model] The value of the logical flag, if the model should
    be scaled or not. (Default FALSE)
\item[adjust.for.con.comp] The value of the logical flag if the
    \texttt{constr=TRUE} option should be reinterpreted.
\end{description}

\subsubsection*{The case \texttt{(scale.model==FALSE \&\&
        adjust.for.con.comp == FALSE)}}

The option \texttt{constr=TRUE} is interpreted as a sum-to-zero
constraint over the whole graph. Singletons are given a uniform
distribution on $(-\infty,\infty)$ before the constraint.

\subsubsection*{The case \texttt{(scale.model==TRUE \&\&
        adjust.for.con.comp == FALSE)}}

The option \texttt{constr=TRUE} is interpreted as a sum-to-zero
constraint over the whole graph. Let $Q = \tau R$ be the standard
precision matrix from the \texttt{besag}-model with precision
parameter $\tau$. Then $R$, except the singletons, are scaled so that
the geometric mean of the marginal variances is 1, and $R$ is modified
so that singletons have a standard Normal distribution.

\subsubsection*{The case \texttt{(scale.model==FALSE \&\&
        adjust.for.con.comp == TRUE)}}

The option \texttt{constr=TRUE} is interpreted as one sum-to-zero
constraint over each of the \texttt{nc2} connected components of size
$\ge2$. Singletons are given a uniform distribution on
$(-\infty,\infty)$.

\subsubsection*{The case \texttt{(scale.model==TRUE \&\&
        adjust.for.con.comp == TRUE)}}

The option \texttt{constr=TRUE} is interpreted as \texttt{nc2}
sum-to-zero constraints for each of the connected components of size
$\ge2$. Let $Q = \tau R$ be the standard precision matrix from the
\texttt{besag}-model with precision parameter $\tau$. Then $R$, are
scaled so that the geometric mean of the marginal variances in each
connected component of size $\ge2$ is 1, and modified so that
singletons have a standard Normal distribution.

%%\end{document}

% LocalWords:  nc constr besag
