#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#define MAX_CITIES 100
#define RADIUS_EARTH 6371
#define PI 3.14159265358979323846

// Deklarasi Struct City untuk mengambil nama kota, latitude, dan longitude
typedef struct {
    char name[50];
    double latitude;
    double longitude;
} City;

// Deklarasi fungsi
double degreesToRadians(double degrees);
double haversine(double lat1, double lon1, double lat2, double lon2);
void swap(int *a, int *b);
double currentRouteDistance(City cities[], int route[], int totalCities);
void bruteForceHelper(City cities[], int route[], int location, int totalCities, double *minDistance, int *bestRoute);
void bruteForceAlgorithm(City cities[], int totalCities, int startLocation);
void DFShelper(double **adj, int n, int start, bool *visited, int *path, int *final_path, double *final_res, double curr_weight, int level);
void depthFirstSearchAlgorithm(City cities[], int totalCities, int startLocation);
void TSPRec(double **adj, double curr_bound, double curr_weight, int level, int curr_path[], int n);
void branchAndBoundAlgorithm(City cities[], int totalCities, int startLocation);
void runAlgorithm(int choice, City cities[], int totalCities, int startLocation);

// Algoritma proses Haversine
// Fungsi degreestoRadians
double degreesToRadians(double degrees) {
    return degrees * PI / 180.0;
}

// Fungsi Haversine
double haversine(double lat1, double lon1, double lat2, double lon2) {
    double dLat = degreesToRadians(lat2 - lat1);
    double dLon = degreesToRadians(lon2 - lon1);
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(degreesToRadians(lat1)) * cos(degreesToRadians(lat2)) *
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = RADIUS_EARTH * c;
    return distance;
}

// Algoritma Brute Force
// Fungsi swap
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Fungsi currentRouteDistance
double currentRouteDistance(City cities[], int route[], int totalCities) {
    double totalDistance = 0.0;
    for (int i = 0; i < totalCities - 1; i++) {
        totalDistance += haversine(cities[route[i]].latitude, cities[route[i]].longitude,
                                   cities[route[i + 1]].latitude, cities[route[i + 1]].longitude);
    }
    totalDistance += haversine(cities[route[totalCities - 1]].latitude, cities[route[totalCities - 1]].longitude,
                               cities[route[0]].latitude, cities[route[0]].longitude);
    return totalDistance;
}

// Fungsi bruteForceHelper
void bruteForceHelper(City cities[], int route[], int location, int totalCities, double *minDistance, int *bestRoute) {
    if (location == totalCities) {
        double totalDistance = currentRouteDistance(cities, route, totalCities);
        if (totalDistance < *minDistance) {
            *minDistance = totalDistance;
            memcpy(bestRoute, route, totalCities * sizeof(int));
        }
        return;
    }

    for (int i = location; i < totalCities; i++) {
        swap(&route[location], &route[i]);
        double partialDistance = currentRouteDistance(cities, route, location + 1);
        if (partialDistance < *minDistance) {  // Use pointer to minDistance
            bruteForceHelper(cities, route, location + 1, totalCities, minDistance, bestRoute);
        }
        swap(&route[location], &route[i]);
    }
}

// Fungsi bruteForceAlgorithm
void bruteForceAlgorithm(City cities[], int totalCities, int startLocation) {
    int route[MAX_CITIES];
    double minDistance = INFINITY;
    int bestRoute[MAX_CITIES];

    for (int i = 0; i < totalCities; i++) {
        route[i] = i;
    }

    swap(&route[0], &route[startLocation]);

    clock_t start = clock();
    bruteForceHelper(cities, route, 1, totalCities, &minDistance, bestRoute); // Pass address of minDistance
    clock_t end = clock();
    double timeElapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Best route found:\n%s -> ", cities[startLocation].name);
    for (int i = 1; i < totalCities; i++) {
        printf("%s -> ", cities[bestRoute[i]].name);
    }
    printf("%s\n", cities[startLocation].name);
    printf("Best route distance: %.6f km\n", minDistance);
    printf("Time elapsed: %.6f seconds\n", timeElapsed);
}

// Algoritma DFS (Depth First Search)
// Fungsi DFShelper
void DFShelper(double **adj, int n, int start, bool *visited, int *path, int *final_path, double *final_res, double curr_weight, int level) {
    visited[start] = true;
    path[level] = start;

    if (level == n - 1) { // Jika semua titik sudah dikunjungi
        if (adj[path[level]][path[0]] != 0) { // Pengecekkan apakah ada edge kembali ke titik awal
            double curr_res = curr_weight + adj[path[level]][path[0]]; // Hitung jarak total

            if (curr_res < *final_res) { // Jika jarak total lebih kecil dari jarak terpendek yang sebelumnya
                *final_res = curr_res; // Update jarak terpendek
                memcpy(final_path, path, (n + 1) * sizeof(int)); // Salin path ke final_path
                final_path[n] = path[0]; // Tambah titik awal di akhir path
            }
        }
    } else {
        // Melakukan proses DFS ke semua titik yang belum dikunjungi
        for (int i = 0; i < n; i++) {
            if (!visited[i] && adj[start][i] != 0) {
                double temp_weight = curr_weight + adj[start][i]; 
                if (temp_weight < *final_res) { // Hanya lanjutkan jika berat sementara kurang dari jarak terpendek yang sudah ditemukan
                    DFShelper(adj, n, i, visited, path, final_path, final_res, temp_weight, level + 1);
                }
            }
        }
    }
    // Reset visited setelah selesai menjelajahi semua jalur dari titik ini
    visited[start] = false;
}

// Fungsi depthFirstSearchAlgorithm
void depthFirstSearchAlgorithm(City cities[], int totalCities, int startLocation) {
    double **matriksjarak = (double **)malloc(totalCities * sizeof(double *));
    for (int i = 0; i < totalCities; i++) {
        matriksjarak[i] = (double *)malloc(totalCities * sizeof(double));
    }

    for (int i = 0; i < totalCities; i++) {
        for (int j = 0; j < totalCities; j++) {
            matriksjarak[i][j] = haversine(cities[i].latitude, cities[i].longitude, cities[j].latitude, cities[j].longitude);
        }
    }

    // Inisialisasi variabel untuk menyimpan jalur dan jarak terpendek
    int *final_path = (int *)malloc((totalCities + 1) * sizeof(int));
    bool *visited = (bool *)malloc(totalCities * sizeof(bool));
    int *path = (int *)malloc((totalCities + 1) * sizeof(int));
    double final_res = DBL_MAX;

    // Set semua elemen visited ke false
    for (int i = 0; i < totalCities; i++) {
        visited[i] = false;
    }

    // Memanggil fungsi DFS helper untuk membantu pencarian jalur terpendek
    clock_t start = clock();
    DFShelper(matriksjarak, totalCities, startLocation, visited, path, final_path, &final_res, 0, 0);
    clock_t end = clock();
    double timeElapsed = (double)(end - start) / CLOCKS_PER_SEC;

    // Mencetak Rute yang ditempuh
    printf("Best route found: \n");
    for (int i = 0; i < totalCities; i++) {
        printf("%s -> ", cities[final_path[i]].name);
    }
    printf("%s", cities[startLocation].name);
    printf("\n");
    printf("Best route distance: %.2f km\n", final_res);
    printf("Time elapsed: %.6f seconds\n", timeElapsed);

    // Free alokasi memori
    free(final_path);
    free(visited);
    free(path);
    for (int i = 0; i < totalCities; i++) {
        free(matriksjarak[i]);
    }
    free(matriksjarak);
}

// Algoritma Branch and Bound
int *final_path; // Simpan array yang berisi path terakhir dari salesman
bool *visited; // Simpan array yang berisi node yang sudah dikunjungi
double final_res = DBL_MAX; // Simpan jarak terpendek dari path yang diambil

// Fungsi copytoFinal
void copyToFinal(int curr_path[], int n) {
    for (int i = 0; i < n; i++)
        final_path[i] = curr_path[i];
    final_path[n] = curr_path[0];
}

// Fungsi firstMin
double firstMin(double **adj, int i, int n) {
    double min = DBL_MAX;
    for (int k = 0; k < n; k++)
        if (adj[i][k] < min && i != k)
            min = adj[i][k];
    return min;
}

// Fungsi secondMin
double secondMin(double **adj, int i, int n) {
    double first = DBL_MAX, second = DBL_MAX;
    for (int j = 0; j < n; j++) {
        if (i == j)
            continue;

        if (adj[i][j] <= first) {
            second = first;
            first = adj[i][j];
        } else if (adj[i][j] <= second && adj[i][j] != first)
            second = adj[i][j];
    }
    return second;
}

// Fungsi TSPrec
void TSPRec(double **adj, double curr_bound, double curr_weight,
            int level, int curr_path[], int n) {
    if (level == n) {
        if (adj[curr_path[level - 1]][curr_path[0]] != 0) {
            double curr_res = curr_weight + adj[curr_path[level - 1]][curr_path[0]];
            if (curr_res < final_res) {
                copyToFinal(curr_path, n);
                final_res = curr_res;
            }
        }
        return;
    }

    for (int i = 0; i < n; i++) {
        if (adj[curr_path[level - 1]][i] != 0 && visited[i] == false) {
            double temp = curr_bound;
            curr_weight += adj[curr_path[level - 1]][i];
            if (level == 1)
                curr_bound -= ((firstMin(adj, curr_path[level - 1], n) +
                                firstMin(adj, i, n)) / 2);
            else
                curr_bound -= ((secondMin(adj, curr_path[level - 1], n) +
                                firstMin(adj, i, n)) / 2);
            if (curr_bound + curr_weight < final_res) {
                curr_path[level] = i;
                visited[i] = true;
                TSPRec(adj, curr_bound, curr_weight, level + 1, curr_path, n);
            }
            curr_weight -= adj[curr_path[level - 1]][i];
            curr_bound = temp;
            memset(visited, false, n * sizeof(bool));
            for (int j = 0; j <= level - 1; j++)
                visited[curr_path[j]] = true;
        }
    }
}

// Fungsi branchAndBoundAlgorithm
void branchAndBoundAlgorithm(City cities[], int totalCities, int startLocation) {
    double **adj = (double **)malloc(totalCities * sizeof(double *));
    for (int i = 0; i < totalCities; i++) {
        adj[i] = (double *)malloc(totalCities * sizeof(double));
    }

    for (int i = 0; i < totalCities; i++) {
        for (int j = 0; j < totalCities; j++) {
            adj[i][j] = haversine(cities[i].latitude, cities[i].longitude, cities[j].latitude, cities[j].longitude);
        }
    }

    int *curr_path = (int *)malloc((totalCities + 1) * sizeof(int));
    final_path = (int *)malloc((totalCities + 1) * sizeof(int));
    visited = (bool *)malloc(totalCities * sizeof(bool));

    double curr_bound = 0;
    memset(curr_path, -1, (totalCities + 1) * sizeof(int));
    memset(visited, 0, totalCities * sizeof(bool));

    for (int i = 0; i < totalCities; i++)
        curr_bound += (firstMin(adj, i, totalCities) + secondMin(adj, i, totalCities));
    curr_bound = (curr_bound / 2.0);

    visited[startLocation] = true;
    curr_path[0] = startLocation;

    clock_t start = clock();
    TSPRec(adj, curr_bound, 0, 1, curr_path, totalCities);
    clock_t end = clock();
    double timeElapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Best route found: \n");
    for (int i = 0; i < totalCities; i++) {
        printf("%s -> ", cities[final_path[i]].name);
    }
    printf("%s", cities[startLocation].name);
    printf("\n");
    printf("Best route distance: %.2f km\n", final_res);
    printf("Time elapsed: %.6f seconds\n", timeElapsed);

    free(curr_path);
    free(visited);
    free(final_path);
    for (int i = 0; i < totalCities; i++) {
        free(adj[i]);
    }
    free(adj);
}

// Run Algoritma berdasarkan pilihan
// Fungsi runAlgorithm
void runAlgorithm(int choice, City cities[], int totalCities, int startLocation) {
    switch (choice) {
        case 1:
            bruteForceAlgorithm(cities, totalCities, startLocation);
            break;
        case 2:
            depthFirstSearchAlgorithm(cities, totalCities, startLocation);
            break;
        case 3:
            branchAndBoundAlgorithm(cities, totalCities, startLocation);
            break;
        default:
            printf("Invalid choice\n");
            break;
    }
}

// Fungsi main
int main() {
    City cities[MAX_CITIES];
    int totalCities = 0;
    char filename[100];
    char startCityName[50];
    int startLocation = -1;
    int choice;

    // Input User file eksternal kota dan Kota awal
    printf("Enter file name containing list of cities: ");
    scanf("%s", filename);
    printf("Enter the starting city: ");
    scanf("%s", startCityName);

    // Membuka file eksternal
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file %s, file not found.\n", filename);
        return 1;
    }

    // Parsing file eksternal
    while (fscanf(file, "%[^,],%lf,%lf\n", cities[totalCities].name, &cities[totalCities].latitude, &cities[totalCities].longitude) != EOF) {
        if (strcmp(cities[totalCities].name, startCityName) == 0) {
            startLocation = totalCities;
        }
        totalCities++;
    }

    fclose(file);

    // Pengecekkan apakah kota terdapat pada struct
    if (startLocation == -1) {
        printf("Starting city '%s' not found in the list.\n", startCityName);
        return 1;
    }

    // Input pilihan algoritma
    printf("Choose the algorithm to run:\n");
    printf("1. Brute Force\n");
    printf("2. Depth-First Search\n");
    printf("3. Branch and Bound\n");
    printf("Enter your choice: ");
    scanf("%d", &choice);

    runAlgorithm(choice, cities, totalCities, startLocation);

    return 0;
}
