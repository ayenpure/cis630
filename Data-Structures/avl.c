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

int height(treenode *node) {
  return 0;
  
}

void insertnode(treenode* node,const char* string) {
  if ( strcmp(node->content, string) <= 0 ) {
    //right
    if(node->right == NULL) {
      node->right = newtreenode(string);
    } else
      insertnode(node->right, string);
  } else {
    //left
    if(node->left == NULL) {
      node->left = newtreenode(string);
    } else
      insertnode(node->left, string);
  }
}

treenode* constructtree(const char** fortree) {
  // Construct root
  treenode *root = newtreenode(fortree[0]);
  for(int i = 1; i < NODE_COUNT; i++) {
      insertnode(root, fortree[i]);
  }
  return root;
}


/*treenode* constructtree(const char** fortree) {
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
}*/

void inorder(treenode *node) {
  if(node->left != NULL) {
    inorder(node->left);
  }
  printf("%s \t", node->content);
  if(node->right != NULL) {
    inorder(node->right);
  }
}

void preorder(treenode *node) {
  printf("%s \t", node->content);
  if(node->left != NULL) {
    preorder(node->left);
  }
  if(node->right != NULL) {
    preorder(node->right);
  }
}

void postorder(treenode *node) {
  if(node->left != NULL) {
    postorder(node->left);
  }
  if(node->right != NULL) {
    postorder(node->right);
  }
  printf("%s \t", node->content);
}


int main(int *argc, char** argv) {
  const char *fortree[100] = {
    "koenigsegg","lamborghini","maserati","pagani","tesla",
    "ferrari","alfa romeo","volkswagen","audi","bmw"
  };
  treenode *root = constructtree(fortree);
  inorder(root);
  printf("\n");
  preorder(root);
  printf("\n");
  postorder(root);
  printf("\n");
  return 0;
}
