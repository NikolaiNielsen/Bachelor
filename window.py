import sys
from PyQt5 import QtWidgets, QtCore


class window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("yoyoy")
        self.home()

    def home(self):
        btn = QtWidgets.QPushButton("Quitter", self)
        btn.clicked.connect(self.close_application)
        btn.resize(btn.minimumSizeHint())
        btn.move(0, 0)
        self.show()

    def close_application(self):
        print("Whoa!")
        sys.exit()


def main():
    app = QtWidgets.QApplication(sys.argv)
    GUI = window()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
