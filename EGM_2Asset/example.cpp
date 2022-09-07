#include <windows.h>
#include <stdio.h>

void CreateFolder(const char *path)
{
    if (!CreateDirectory(path, NULL))
    {
        return;
    }
}

int main()
{
    CreateFolder(".\\result\\");
}
