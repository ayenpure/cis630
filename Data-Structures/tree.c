#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define NODE_COUNT 10

typedef struct treenode {
  char *content;
  struct treenode *left;
  struct treenode *right;
} treenode;

treenode* newtreenode(const char* content) {
  treenode* newnode = (treenode*)malloc(sizeof(treenode));
  newnode->content = (char*)malloc(100*sizeof(char));
  strcpy(newnode->content, content);
  newnode->left = NULL;
  newnode->right = NULL;
  return newnode;
}

treenode* constructtree(const char** fortree) {
  // Construct root
  treenode *root = newtreenode(fortree[0]);
  for(int i = 1; i < NODE_COUNT; i++) {
    treenode *currnode = root;
    int none = 1, greater = 0;
    while(none) {
      if(strcmp(currnode->content, fortree[i]) <= 0) {
        if(currnode->right != NULL)
          currnode = currnode->right;
        else {
          none = 0;
          greater = 1;
          break;
        }
      } else {
        if(currnode->left != NULL)
          currnode = currnode->left;
        else {
          none = 0;
          greater = 0;
          break;
        }
      }
    }
    treenode* newnode = newtreenode(fortree[i]);
    if(greater) {
      currnode->right = newnode;
    } else {
      currnode->left = newnode;
    }
  }
  return root;
}

int main(int *argc, char** argv) {
  const char *fortree[100] = {
    "koenigsegg","lamborghini","maserati","pagani","tesla",
    "ferrari","alfa romeo","volkswagen","audi","bmw"
  };
  treenode *content = constructtree(fortree);
  return 0;
}
