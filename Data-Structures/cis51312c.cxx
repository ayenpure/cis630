#include <stdio.h>

#define MAX_NODES 7

typedef struct vertex {
  int label;
  double bottleneck;
}vertex;

int isempty(vertex *vertices) {
  int empty = 1;
  int max = -1;
  int index = -1;
  for(int i = 0;i < MAX_NODES; i++) {
    if(vertices[i].label == -1)
      continue;
    else {
        empty = 0;
        if(vertices[i].bottleneck > max) {
          index = i;
          max = vertices[i].bottleneck;
        }
    }
  }
  if(!empty)
    return index;
}

int main(int argc, char **argv) {
  double G[MAX_NODES][MAX_NODES] = {
    {-1,2,3,4,-1,-1,-1},
    {-1,-1,-1,-1,5,-1,-1},
    {-1,-1,-1,-1,6,5,-1},
    {-1,-1,-1,-1,-1,7,-1},
    {-1,-1,-1,-1,-1,-1,8},
    {-1,-1,-1,-1,-1,-1,9},
    {-1,-1,-1,-1,-1,-1,-1}
  };

  vertex vertices[MAX_NODES];
  for(int i=0;i<MAX_NODES;i++) {
    vertices[i].label = i;
    vertices[i].bottleneck = -1;
  }
  vertices[0].bottleneck = 1000000;
  int index = -1;
  while((index = isempty(vertices)) != -1) {
      vertices[index].label = -1;
      for(int i = 0; i < MAX_NODES; i++) {
        if(G[index][i] != -1) {
          double capacity = G[index][i];
          if (vertices[index].bottleneck > vertices[i].bottleneck) {
            if(capacity <= vertices[index].bottleneck) {
              vertices[i].bottleneck = capacity;
            } else {
              vertices[i].bottleneck = vertices[index].bottleneck;
            }
          }
        }
      }
  }
  return 0;
}

/*if(vertices[i].bottleneck > vertices[i].bottleneck)
continue;
else if(vertices[i].bottleneck > capacity)
continue;
else {
vertices[i].bottleneck = capacity;
}*/


/*if(G[index][i] != -1) {
  double capacity = G[index][i];
  if(vertices[i].bottleneck >= vertices[index].bottleneck){
    continue;
  } else if (vertices[index].bottleneck > vertices[i].bottleneck) {
    if(capacity <= vertices[i].bottleneck)
      continue;
    else if(capacity <= vertices[index].bottleneck) {
      vertices[i].bottleneck = capacity;
    } else {
      vertices[i].bottleneck = vertices[index].bottleneck;
    }
  }
}*/
