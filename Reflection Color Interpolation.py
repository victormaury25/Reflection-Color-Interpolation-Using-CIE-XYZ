# -*- coding: utf-8 -*-
"""
@author: victor
"""

import os
import colour
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay

def XYZ_to_xyY(L):
    """
    converted color from the XYZ space to the xyY space

    Parameters
    ----------
    L : 3D array
        color in the XYZ space

    Returns
    -------
    3D array
        color in the xyY space

    """
    if (L[0] + L[1] + L[2]) == 0:
        return 0, 0, L[1]  # We don't divide by 0

    x = L[0] / (L[0] + L[1] + L[2])
    y = L[1] / (L[0] + L[1] + L[2])

    return [x, y, L[1]]

def xyY_to_XYZ(L):
    """
    converted color from the xyY space to the XYZ space

    Parameters
    ----------
    L : 3D array
        color in the xyY space

    Returns
    -------
    3D array
        color in the XYZ space

    """
    if L[1] == 0:
        return 0, 0, 0  # we don't divide by 0

    X = (L[0] / L[1]) * L[2]
    Z = ((1 - L[0] - L[1]) / L[1]) * L[2]

    return [X, L[2], Z]


def xyz_to_rgb(XYZ):
    """
    Convert XYZ values to RGB values in the sRGB color space.

    Parameters:
        XYZ :array
            XYZ values.

    Returns:
        RGB : array
            Gamma-corrected RGB values.
    """
    M = np.array([
        [3.2406, -1.5372, -0.4986],
        [-0.9689, 1.8758, 0.0415],
        [0.0557, -0.2040, 1.0570]
    ])

    RGB = np.dot(M, XYZ)
    RGB = np.clip(RGB, 0, 1)  # Clip values to [0, 1]

    # Apply gamma correction
    gamma_corrected = np.where(RGB <= 0.0031308,
                               12.92 * RGB,
                               1.055 * (RGB ** (1 / 2.4)) - 0.055)
    return gamma_corrected


# Function to convert spectral reflectance values to CIE XYZ values
def passage_spectre_XYZ(valeurs_alb, angle_source, angle_capteur, n, ks):
    """
    Convert spectral reflectance values to XYZ values using the CIE 1931 2 Degree Standard Observer.

    Parameters:
        valeurs_alb : array
            Reflectance values corresponding to wavelengths.
        angle_source : float
            angle of the light source with the normal of the reflexion surface in degrees
        angle_capteur : float
            angle of the captor with the normal of the reflexion surface in degrees
        n : int
            power of the cosinus we want to use for our Phong BRDF model
        ks : float
            value of the specular coefficient

    Returns:
        XYZ (array): Calculated XYZ values.
    """
    cmfs = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
    wavelengths = cmfs.wavelengths
    mask = (wavelengths >= 380) & (wavelengths <= 749)  # Filter wavelengths from 380 nm to 749 nm
    wavelengths = wavelengths[mask]
    k_s=np.array([ks]*len(valeurs_alb))
    valeurs_alb=np.array(valeurs_alb)
    albedo=valeurs_alb+k_s*np.cos((angle_capteur-angle_source)*np.pi/180)**n

    if len(wavelengths) == len(albedo):
        X_sensitivity = cmfs.values[mask, 0]
        Y_sensitivity = cmfs.values[mask, 1]
        Z_sensitivity = cmfs.values[mask, 2]

        # Compute XYZ values using weighted integration
        X = np.trapz(albedo * X_sensitivity, wavelengths)
        Y = np.trapz(albedo * Y_sensitivity, wavelengths)
        Z = np.trapz(albedo * Z_sensitivity, wavelengths)

        XYZ = np.array([X, Y, Z]) * 0.0113  # Normalize XYZ values
        return XYZ
    else:
        raise ValueError(
            f"Length mismatch: wavelengths ({len(wavelengths)}) and valeurs_alb ({len(valeurs_alb)})."
        )

# Function to read albedo values from a CSV file
def lire_albedo_csv(csv_file):
    """
    Read albedo values from a CSV file.

    Parameters:
        csv_file (str): Path to the CSV file.

    Returns:
        albedo_values (array): Albedo values from the second column.
    """
    data = pd.read_csv(csv_file)
    albedo_values = data.iloc[:, 1].values  # Assume the second column contains albedo values
    return albedo_values

liste_points_sources_XYZ=np.array([[1.939,1,0.091],
                          [0.5,1,0.167],
                          [2.5,1,13.167],
                          [0.880,1,0.320],
                          [0.667,1,2],
                          [1.667,1,5.333],
                          [1.250,1,0.750],
                          [1,1,1.030],
                          [0.286,1,2.457],
                          [1.833,1,0.583]]) #values of each vertix of the mesh in the XYZ space of color


def interpoler_delaunay(point, liste_pts_sources, valeurs_associees):
    """
    Interpolation in 2 dimensions in a Delaunay mesh

    Parameters
    ----------
    point : 2D array
        coorodinates of the point we want to interpolate.
    liste_pts_sources : N-D array
        list of all the points used for our Delaunay mesh (here each values are 2D arrays representing the xy coordinates)
    valeurs_associees : N-D array
        list of all the values associated to each points of the mesh (here each values are 2D arrays representing the xy coordinates)

    Returns
    -------
    2D array
        interpolated values and the weight of each vertices of the triangle our point is inside of

    """
    triangulation = Delaunay(liste_pts_sources)
    # Trouver le simplex contenant le point
    simplex_index = triangulation.find_simplex(point)
    if simplex_index == -1:
        return None, None  # Le point est en dehors de la triangulation

    # Obtenir les indices des sommets du triangle
    triangle_vertices = triangulation.simplices[simplex_index]
    # Obtenir les coordonnées des sommets
    vertices_coords = liste_pts_sources[triangle_vertices]

    # Calcul des coordonnées barycentriques
    T = np.hstack((vertices_coords, np.ones((3, 1))))  # Ajout de la coordonnée pour homogénéisation
    pt_h = np.append(point, 1)  # Homogénéisation du point
    bary_coords = np.linalg.solve(T.T, pt_h)

    # Vérifier si le point est bien dans le triangle (en pratique, pas nécessaire avec find_simplex)
    if np.any(bary_coords < 0):
        return None, None

    # Interpolation des valeurs associées (x', y')
    interpolated_values = np.dot(bary_coords, valeurs_associees[triangle_vertices])
    return [interpolated_values, bary_coords]


def display_color(XYZ):
    """
    Display a color patch for given XYZ values.

    Parameters:
        XYZ (array): XYZ values.
    """
    RGB = xyz_to_rgb(XYZ)
    plt.figure(figsize=(2, 2))
    plt.gca().set_facecolor("black")
    plt.imshow([[RGB]], extent=(0, 1, 0, 1))
    plt.title(f"Color (XYZ = {np.round(XYZ, 3)})")
    plt.axis("off")
    plt.show()


def couleur_reflexion(albedo_csv, couleur_source_XYZ, angle_source, angle_capteur, n, ks):
    """
    Display the color and XYZ values of the relexion on a given material exposed to a given light source

    Parameters
    ----------
    albedo_csv : string
        the path to the csv file with the albedo values for the material studied
    couleur_source_XYZ : array
        XYZ coordinates of the color of the light source
    angle_source : float
        angle of the light source with the normal of the reflexion surface in degrees
    angle_capteur : float
        angle of the captor with the normal of the reflexion surface in degrees
    n : int
        power of the cosinus we want to use for our Phong BRDF model
    ks : float
        value of the specular coefficient
    """
    albedo = lire_albedo_csv(albedo_csv)
    albedo_XYZ = passage_spectre_XYZ(albedo, angle_source, angle_capteur, n, ks)
    
    D=[]
    L=[]
    for i in range(len(liste_points_sources_XYZ)):
        D.append(albedo_XYZ*liste_points_sources_XYZ[i])
        L.append([XYZ_to_xyY(liste_points_sources_XYZ[i])[0],XYZ_to_xyY(liste_points_sources_XYZ[i])[1]])

    liste_valeurs_associees_XYZ=np.array(D)#values associated to each vertix of the mesh converted in the XYZ space of color
    liste_points_sources_xy=np.array(L)#xy values of each vertix of the mesh converted in the xy space of color

    M=[]
    for i in range(len(liste_valeurs_associees_XYZ)):
        M.append([XYZ_to_xyY(liste_valeurs_associees_XYZ[i])[0],XYZ_to_xyY(liste_valeurs_associees_XYZ[i])[1]])

    liste_valeurs_associees_xy=np.array(M) #values associated to each vertix of the mesh converted in the xy space of color

    point_XYZ=couleur_source_XYZ
    point_xy=[XYZ_to_xyY(point_XYZ)[0],XYZ_to_xyY(point_XYZ)[1]]

    couleur_reflexion_xy= interpoler_delaunay(point_xy, liste_points_sources_xy, liste_valeurs_associees_xy)[0]
    couleur_reflexion_xyY=np.append(couleur_reflexion_xy, albedo_XYZ[1]*point_XYZ[1])
    couleur_reflexion_XYZ=xyY_to_XYZ(couleur_reflexion_xyY)
    display_color(couleur_reflexion_XYZ) 