class Classname:
    @staticmethod
    def fun():
        print('静态方法')

    @classmethod
    def a(cls):
        print('类方法')

    # 普通方法
    def b(self):
        print('普通方法')



Classname.fun()
Classname.a()

C = Classname()
C.fun()
C.a()
C.b()
