from __future__ import print_function
import sys
import argparse
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import newton
import numpy as np

def estimate_threshold(network_size, covariance):
    """Estimate the threshold based upon the configuration model."""
    func = lambda t: lambda y: (1-erf((t-y*covariance)/
        (np.sqrt(2*(1-covariance**2)))))*np.exp(-y**2/2)
    func2 = lambda y: np.exp(-y**2/2)
    q2 =  lambda t: (
        2*np.sqrt(2*np.pi))**(-1)*quad(func(t), t, np.inf)[0]
    q1 = lambda t: (1/np.sqrt(2*np.pi))*quad(func2, t, np.inf)[0]
    func3 =  lambda t: network_size*q2(t)/(q1(t))-1
    return newton(func3, 0.1)

def get_mean_degree(network_size, threshold):
    """Get the expected degree for a node."""
    func = lambda y: np.exp(-y**2/2)
    q1 = lambda t: (1/np.sqrt(2*np.pi))*quad(func, t, np.inf)[0]
    return (network_size-1)*q1(threshold)

def main(arguments):
    """Estimate the critical threshold for the giant component
    transition, for a fixed covariance. Then return a list of values
    around it to perform simulations.
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-N', '--network_size', type=int,
                        help="Number of nodes")
    parser.add_argument('-c', '--covariance', type=float,
                        help="Correlation parameter")
    parser.add_argument('-n', '--list_size', type=int, help="""Size of the list
                        of values around the critical threshold""")
    parser.add_argument('-f', '--fudge_factor', type=float, default=1,
                        help="Factor to adjust the range due to discrepancy")
    args = parser.parse_args(arguments)

    #determine the critical threshold estimated by the CM
    network_size =  args.network_size
    covariance = args.covariance
    critical_threshold = estimate_threshold(network_size,covariance)

    #fudge the threshold estimate due to clustering effects
    null_covariance_threshold = estimate_threshold(network_size, 0.)
    threshold_estimate = (null_covariance_threshold +
                          (1-covariance)*args.fudge_factor*
                          (critical_threshold - null_covariance_threshold))

    #create a list of values around it
    threshold_list = np.linspace(threshold_estimate - 0.5,
                                 threshold_estimate + 0.5,
                                 args.list_size)
    threshold_list_string = ""
    for t in threshold_list:
        threshold_list_string += "{0:.5f} ".format(t)

    #print the resulting range of values in a string
    print(threshold_list_string[:-1])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
