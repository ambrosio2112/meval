
#' @title Modelo de evaluaciones para individuos basado en Val Jonhson
#'
#' @description
#' Modelo de evaluaciones para individuos basado en Val Jonhson.
#' Este modelo estima por MCMC los coeficientes de la relación entre las evaluaciones y y las covariables X
#' El modelo soporta que los individuos no evaluen todos las mismas marcas (por ejemplo cuando no las conocen)
#' Adicionalmente, el modelo puede estimar los coeficientes en nclases grupos. Si nclases=N el numero de individuos, equivale el modelo jerarquico por individuos,
#' pero nclases puede ser menor, por ejemplo, coeficientes por estrato.
#'
#' @param N Numero de individuos
#' @param J Numero de objetos (marcas, atributos) evaluados
#' @param E Valor maximo de la escala de evaluacione 1,..,E
#' @param K Cantidad de coeficientes a estimar (nbeta)
#' @param nclases Cantidad de clases a estimar
#' @param clases  Mapa entre las filas de individuos y las clases
#' @param marcaevaluada Indicador NxJ de marcas evaluadas por individuo
#' @param X Covariables NxJxK
#' @param y Evaluaciones realizadas NxJ
#' @param niter Numero de iteraciones
#' @param nburn Cantidad de iteraciones a quemar del MCMC
#' @param nstep Pasos para tomar la simulacion del MCMC
#' @param nadjust No recuerdo
#' @param alpha Parametros de las previas
#' @param lambda Parametros de las previas
#' @param tau0 Parametros de las previas
#' @param sigmaMH Parametros de las previas
#'
#' @return El programa genera archivos de texto con las simulaciones realizadas
#' @export
#' @useDynLib eval
#' @examples
#' ## Falta buscar el ejemplo, por ejemplo en los proyectos de modelos de marcas o jerarquizacion de drivers.
#'
evaluaciones <- function(N,J,E,K,nclases, clases,
                         marcaevaluada,X,y,
                         niter=15000, nburn=5000, nstep=100, nadjust=100,
                         alpha = as.double(0.5),
                         lambda = as.double(2),
                         tau0 = as.double(2),
                         sigmaMH = as.double(1))
  {
    z <- .Fortran("evaluaciones",
                  niter=as.integer(niter),
                  nburn=as.integer(nburn),
                  nstep=as.integer(nstep),
                  nadjust=as.integer(nadjust),
                  N=as.integer(N),
                  nmarcas=as.integer(J),
                  E=as.integer(E),
                  nbeta=as.integer(K),
                  nclass=as.integer(nclases),
                  marcaevaluada=as.integer(marcaevaluada),
                  ev=as.double(X),
                  y=as.integer(y),
                  missing=as.integer(missing),
                  class=as.integer(clases),
                  alpha=as.double(alpha),
                  lambda=as.double(lambda),
                  sigmaMH=as.double(sigmaMH),
                  tau0=as.double(tau0),
                  Sbeta=as.double(rep(2*tau0,K)),
                  mubeta=as.double(rep(0,K)),
                  NAOK=TRUE
                  )
    NULL
  }
